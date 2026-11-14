"""PPI + target-gene evaluation pipeline on synthetic fixtures (no biology implied).

UniProt and the PLM-interact model are replaced by deterministic fakes; the code
paths around them (validation, provenance, refusal rules) are the real ones.
"""
import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

import analyze
import commot_baseline
import evaluate_targets as ev
import plm_score
import proteins
import research
import stage1

LETTERS = 'ACDEFGHIKLMNPQRSTVWY'


def quiet():
    return contextlib.redirect_stdout(io.StringIO())


def feature(kind, start, end, description=''):
    return {'type': kind, 'description': description,
            'location': {'start': {'value': start}, 'end': {'value': end}}}


def fake_entry(symbol, accession, length=60, taxon=10090, reviewed=True, seed=0, features=()):
    rng = np.random.default_rng(seed)
    return {'primaryAccession': accession, 'uniProtkbId': f'{symbol.upper()}_MOUSE', 'features': list(features),
            'entryType': 'UniProtKB reviewed (Swiss-Prot)' if reviewed else 'UniProtKB unreviewed (TrEMBL)',
            'organism': {'taxonId': taxon}, 'genes': [{'geneName': {'value': symbol}}],
            'proteinDescription': {'recommendedName': {'fullName': {'value': f'{symbol} protein'}}},
            'entryAudit': {'entryVersion': 7, 'sequenceVersion': 1},
            'sequence': {'value': ''.join(rng.choice(list(LETTERS), length))}}


def fake_uniprot(entries):
    """entries: symbol -> list of entries returned by the search."""
    by_accession = {e['primaryAccession']: e for es in entries.values() for e in es}

    def fetch(url):
        headers = {'x-uniprot-release': '2026_04', 'x-uniprot-release-date': '01-Oct-2026'}
        if '/search?' in url:
            symbol = url.split('gene_exact%3A')[1].split('%29')[0]
            return {'results': entries.get(symbol, [])}, headers
        return by_accession[url.rsplit('/', 1)[1][:-5]], headers
    return fetch


class Pipeline(unittest.TestCase):
    def build(self, root):
        rng = np.random.default_rng(3)
        grid = np.array([[i, j] for i in range(6) for j in range(5)], float)
        xy = grid + rng.uniform(-.2, .2, grid.shape)
        n = len(xy)
        receptor_level = rng.poisson(6, n) + 1
        genes = ['Cxcl12', 'Cxcr4', 'Cxcr3', 'Dll1', 'Notch1', 'Notch2'] + [f'Tg{i}' for i in range(6)] + ['Other']
        x = np.column_stack([rng.poisson(8, n) + 1, receptor_level, rng.poisson(5, n) + 1,
                             rng.poisson(6, n) + 1, rng.poisson(5, n) + 1, rng.poisson(4, n) + 1]
                            + [rng.poisson(receptor_level + 2) for _ in range(6)] + [rng.poisson(50, n)])
        data = ad.AnnData(sparse.csr_matrix(x.astype(float)), obs=pd.DataFrame(index=[f's{i}' for i in range(n)]),
                          var=pd.DataFrame(index=genes))
        data.obsm['spatial'] = xy
        data.obs['scc_anno'] = rng.choice(['A', 'B'], n)
        data.write_h5ad(root/'fixture.h5ad')
        definitions = {'diff_LigRecDB.csv': ('diffusion', [('CXCL12', 'CXCR4'), ('CXCL12', 'CXCR3')]),
                       'cont_LigRecDB.csv': ('contact', [('DLL1', 'NOTCH1'), ('DLL1', 'NOTCH2')])}
        (root/'lr').mkdir()
        checksums = {}
        for filename, (kind, rows) in definitions.items():
            path = root/'lr'/filename
            pd.DataFrame(rows, columns=['source', 'target']).assign(origin='fixture').to_csv(path, index=False)
            checksums[filename] = (kind, stage1.digest(path))
        args = stage1.parser().parse_args(['prepare', '--h5ad', str(root/'fixture.h5ad'), '--lr-dir', str(root/'lr'),
            '--output', str(root/'input'), '--matrix', 'X', '--expression-state', 'counts',
            '--symbol-policy', 'author-uppercase', '--coordinate-unit', 'synthetic', '--max-spots', str(n),
            '--pairs-per-type', '2', '--min-fraction', '0', '--min-mean', '0'])
        return args, checksums

    def test_end_to_end(self):
        with tempfile.TemporaryDirectory() as d, quiet():
            root = Path(d)
            args, checksums = self.build(root)
            with patch.dict(stage1.LR_FILES, checksums, clear=True):
                stage1.prepare(args)
                symbols = ['Cxcl12', 'Cxcr4', 'Cxcr3', 'Dll1', 'Notch1', 'Notch2']
                fetch = fake_uniprot({s: [fake_entry(s, f'P{k:05d}', 40 + 10*k, seed=k)] for k, s in enumerate(symbols)})
                margs = proteins.parser().parse_args(['map', '--input', str(root/'input'), '--output', str(root/'prot')])
                proteins.map_proteins(margs, fetch=fetch)
                pmap = pd.read_csv(root/'prot'/'protein_map.csv')
                self.assertEqual(dict(zip(pmap.lr_symbol, pmap.data_symbol))['CXCR4'], 'Cxcr4')
                proteins.make_pairs(proteins.parser().parse_args(['pairs', '--input', str(root/'input'),
                    '--proteins', str(root/'prot'), '--output', str(root/'jobs'), '--lr-dir', str(root/'lr')]))
            jobs = pd.read_csv(root/'jobs'/'jobs.csv')
            self.assertEqual((jobs.set == 'cognate').sum(), 8)        # 4 pairs x 2 orders
            self.assertEqual((jobs.set == 'panel_cross').sum(), 8)    # CXCL12xNOTCH*, DLL1xCXCR*
            # Fake model: deterministic score from the sequences; real score() path otherwise.
            checkpoint = root/'pytorch_model.bin'; checkpoint.write_bytes(b'fake weights')
            fake = lambda q, t: ((len(q) * 7 + len(t) * 3) % 97 / 100 + .01, len(q) + len(t) + 3)
            sargs = plm_score.parser().parse_args(['score', '--jobs', str(root/'jobs'/'plm_input.csv'),
                '--checkpoint', str(checkpoint), '--base-model', 'fake/esm', '--checkpoint-name', 'fake/plm',
                '--checkpoint-revision', 'abc', '--output', str(root/'scores.csv')])
            with patch.object(plm_score, 'build_scorer', return_value=(fake, dict(device='cpu'))):
                plm_score.score(sargs)
            proteins.import_scores(proteins.parser().parse_args(['import', '--input', str(root/'input'),
                '--jobs', str(root/'jobs'), '--scores', str(root/'scores.csv'), '--output', str(root/'ppi')]))
            ppi = pd.read_csv(root/'ppi'/'ppi.csv')
            np.testing.assert_allclose(ppi.score, (ppi.score_LR + ppi.score_RL)/2)
            report = json.loads((root/'ppi'/'score_report.json').read_text())
            self.assertIn('fake/plm@abc', report['score_provenance'])
            with patch.dict(stage1.LR_FILES, checksums, clear=True):
                stage1.attach_ppi(stage1.parser().parse_args(['attach-ppi', '--input', str(root/'input'),
                    '--ppi-dir', str(root/'ppi'), '--output', str(root/'input_ppi')]))
                for name in ['expression.csv', 'spatial.csv', 'lr.csv']:
                    self.assertEqual(stage1.digest(root/'input'/name), stage1.digest(root/'input_ppi'/name))
                stage1.run_stage(stage1.parser().parse_args(['run', '--input', str(root/'input'),
                                                             '--output', str(root/'res_s1')]))
                stage1.run_stage(stage1.parser().parse_args(['run', '--input', str(root/'input_ppi'),
                                                             '--output', str(root/'res_ppi')]))
                stage1.export_genes(stage1.parser().parse_args(['export-genes', '--input', str(root/'input'),
                                                                '--output', str(root/'genes')]))
            summary = pd.read_csv(root/'res_ppi'/'ot_00'/'results'/'summary.csv')
            self.assertEqual(set(summary.family), {'r0_raw'})
            self.assertIn('ppi_constant', set(summary.score_mode))
            self.assertTrue(summary.score_mode.str.startswith('perm_diffusion').any())
            controls = pd.read_csv(root/'res_ppi'/'ot_00'/'analysis'/'control_matched.csv')
            self.assertEqual(len(controls), len(summary) - 1)
            self.assertTrue(pd.read_csv(root/'res_ppi'/'ot_00'/'analysis'/'shared_entity_split.csv').size)
            commot_baseline.run(root/'input', root/'commot', dis_thr_h=(3, 5))
            db = root/'db'; db.mkdir()
            pd.DataFrame(dict(receptor=['Cxcr4', 'Cxcr4', 'Cxcr3', 'Notch1', 'Notch2'], pathway='p',
                              tf=['Tf1', 'Tf2', 'Tf2', 'Tf3', 'Tf3'], weight=[.9, .1, .5, .4, .4],
                              source='fixture')).to_csv(db/'receptor_tf.csv', index=False)
            pd.DataFrame(dict(tf=['Tf1']*5 + ['Tf2']*4 + ['Tf3']*6 + ['Tf3'],
                              target=[f'Tg{i}' for i in range(5)] + [f'Tg{i}' for i in range(1, 5)]
                              + [f'Tg{i}' for i in range(6)] + ['Cxcr4'])).to_csv(db/'tf_tg.csv', index=False)
            (db/'db.json').write_text(json.dumps(dict(status='complete', converted_sha256={
                f: ev.digest(db/f) for f in ['receptor_tf.csv', 'tf_tg.csv']})))
            eargs = ev.parser().parse_args(['run', '--input', str(root/'input'), '--genes', str(root/'genes'),
                '--db', str(db), '--output', str(root/'eval'), '--top-tf', '1', '2', '--covariate', 'scc_anno',
                '--results', str(root/'res_s1'/'ot_00'/'results'), str(root/'res_ppi'/'ot_00'/'results'),
                str(root/'commot')])
            ev.evaluate(eargs)
            sets = pd.read_csv(root/'eval'/'target_sets.csv')
            top1 = sets[(sets.receptor == 'CXCR4') & (sets.top_tf.astype(str) == '1')].iloc[0]
            self.assertEqual(top1.tfs, 'Tf1')                       # highest weight only
            self.assertNotIn('Cxcr4', str(sets.targets.tolist()))   # panel genes removed
            metrics = pd.read_csv(root/'eval'/'case_metrics.csv')
            self.assertEqual(set(metrics.family), {'inverse', 'r0_raw', 'commot_euc'})
            self.assertTrue(metrics.spearman.between(-1, 1).all())
            self.assertGreater(len(pd.read_csv(root/'eval'/'permutation_summary.csv')), 0)
            self.assertGreater((root/'eval'/'target_eval.png').stat().st_size, 1000)
            # Results computed on a different LR panel must be refused.
            other = root/'other'; other.mkdir()
            for name in ['expression.csv', 'spatial.csv', 'cell_meta.csv']:
                (other/name).write_bytes((root/'input'/name).read_bytes())
            pd.read_csv(root/'input'/'lr.csv').head(3).to_csv(other/'lr.csv', index=False)
            eargs.input, eargs.output = str(other), str(root/'eval2')
            with self.assertRaisesRegex(ValueError, 'differs from --input'):
                ev.evaluate(eargs)


class Units(unittest.TestCase):
    def test_uniprot_rules(self):
        two = fake_uniprot({'Gene': [fake_entry('Gene', 'P1'), fake_entry('Gene', 'P2')]})
        with self.assertRaisesRegex(ValueError, '2 reviewed mouse entries'):
            proteins.resolve('Gene', fetch=two)
        synonym_only = fake_uniprot({'Gene': [fake_entry('Other', 'P3')]})
        with self.assertRaisesRegex(ValueError, '0 reviewed'):
            proteins.resolve('Gene', fetch=synonym_only)
        human = fake_uniprot({'Gene': [fake_entry('Gene', 'Q9', taxon=9606)]})
        with self.assertRaisesRegex(ValueError, 'not a mouse entry'):
            proteins.resolve('Gene', override='Q9', fetch=human)
        self.assertEqual(proteins.auroc([.9, .8], [.1, .8]), .875)

    def test_map_reports_all_ambiguous_genes_then_accepts_overrides(self):
        with tempfile.TemporaryDirectory() as d, quiet():
            root = Path(d)
            (root/'input').mkdir()
            pd.DataFrame(dict(ligand=['NRXN1', 'CXCL12'], receptor=['NLGN1', 'CXCR4'],
                              type=['contact', 'diffusion'])).to_csv(root/'input'/'lr.csv', index=False)
            table = root/'symbols.csv'
            pd.DataFrame(dict(lr_symbol=['NRXN1', 'NLGN1', 'CXCL12', 'CXCR4'],
                              data_symbol=['Nrxn1', 'Nlgn1', 'Cxcl12', 'Cxcr4'])).to_csv(table, index=False)
            fetch = fake_uniprot({'Nrxn1': [fake_entry('Nrxn1', 'Q9CS84', 300, seed=1),
                                            fake_entry('Nrxn1', 'P0DI97', 90, seed=2)],
                                  'Nlgn1': [fake_entry('Nlgn1', 'Q99K10', 120, seed=3)],
                                  'Cxcl12': [], 'Cxcr4': [fake_entry('Cxcr4', 'P70658', 80, seed=4)]})
            args = proteins.parser().parse_args(['map', '--input', str(root/'input'), '--output',
                                                 str(root/'prot'), '--symbol-table', str(table)])
            with self.assertRaisesRegex(ValueError, r"2 gene\(s\).*\['CXCL12', 'NRXN1'\]"):
                proteins.map_proteins(args, fetch=fetch)
            self.assertFalse((root/'prot').exists())
            candidates = pd.read_csv(root/'prot_candidates.csv')
            self.assertEqual(set(candidates.accession.dropna()), {'Q9CS84', 'P0DI97'})
            self.assertEqual(candidates.loc[candidates.lr_symbol == 'CXCL12', 'n_hits'].item(), 0)
            fetch = fake_uniprot({'Nrxn1': [fake_entry('Nrxn1', 'Q9CS84', 300, seed=1),
                                            fake_entry('Nrxn1', 'P0DI97', 90, seed=2)],
                                  'Nlgn1': [fake_entry('Nlgn1', 'Q99K10', 120, seed=3)],
                                  'Cxcl12': [fake_entry('Cxcl12', 'P40224', 90, seed=5)],
                                  'Cxcr4': [fake_entry('Cxcr4', 'P70658', 80, seed=4)]})
            pd.DataFrame(dict(lr_symbol=['NRXN1', 'CXCL12'], accession=['P0DI97', 'P40224'],
                              reason=['beta form; contains LNS6', 'checked manually'])).to_csv(root/'o.csv', index=False)
            args.override = str(root/'o.csv')
            proteins.map_proteins(args, fetch=fetch)
            pmap = pd.read_csv(root/'prot'/'protein_map.csv').set_index('lr_symbol')
            self.assertEqual(pmap.loc['NRXN1', 'accession'], 'P0DI97')
            self.assertEqual(pmap.loc['NRXN1', 'resolution'], 'manual_override')
            pd.DataFrame(dict(lr_symbol=['Nrxn1'], accession=['P0DI97'], reason=['x'])).to_csv(root/'bad.csv', index=False)
            args.override, args.output = str(root/'bad.csv'), str(root/'prot2')
            with self.assertRaisesRegex(ValueError, 'not in this LR panel'):
                proteins.map_proteins(args, fetch=fetch)
            for region in ['full', 'extracellular']:
                proteins.make_pairs(proteins.parser().parse_args(['pairs', '--input', str(root/'input'),
                    '--proteins', str(root/'prot'), '--output', str(root/f'jobs_{region}'), '--region', region]))
            full = pd.read_csv(root/'jobs_full'/'jobs.csv')
            self.assertTrue(full.query_region.str.startswith('full').all())
            self.assertEqual(json.loads((root/'jobs_extracellular'/'jobs.json').read_text())['region'], 'extracellular')

    def test_extracellular_region_rule(self):
        seq_entry = lambda *f: fake_entry('G', 'P1', 100, features=f)
        type1 = seq_entry(feature('Signal', 1, 20), feature('Topological domain', 21, 70, 'Extracellular'),
                          feature('Transmembrane', 71, 91), feature('Topological domain', 92, 100, 'Cytoplasmic'))
        seq, used = proteins.region_sequence(type1, 'extracellular')
        self.assertEqual((len(seq), used), (50, 'extracellular 21-70'))
        self.assertEqual(seq, type1['sequence']['value'][20:70])
        gpi = seq_entry(feature('Signal', 1, 18), feature('Chain', 19, 80), feature('Propeptide', 81, 100))
        self.assertEqual(proteins.region_sequence(gpi, 'extracellular')[1], 'mature chain 19-80')
        bare = seq_entry(feature('Signal', 1, 25))
        self.assertEqual(proteins.region_sequence(bare, 'extracellular')[1], 'without signal peptide 26-100')
        self.assertEqual(proteins.region_sequence(bare, 'full')[1], 'full 1-100')
        multipass = seq_entry(feature('Topological domain', 1, 10, 'Extracellular'),
                              feature('Topological domain', 40, 50, 'Extracellular'))
        with self.assertRaisesRegex(ValueError, 'multi-pass'):
            proteins.region_sequence(multipass, 'extracellular')

    def test_plm_refuses_long_altered_or_invalid_input(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'jobs.csv'
            pd.DataFrame(dict(pair_id=['a', 'b'], query=['ACD', 'A'*1000], text=['KLM', 'C'*700])).to_csv(path, index=False)
            jobs = plm_score.load_jobs(path)
            good = lambda q, t: (.5, len(q) + len(t) + 3)
            with self.assertRaisesRegex(ValueError, 'exceed 1603'):
                plm_score.run_jobs(jobs, good)
            self.assertEqual(len(plm_score.run_jobs(jobs, good, allow_long=True)[0]), 2)
            with self.assertRaisesRegex(ValueError, 'truncated'):
                plm_score.run_jobs(jobs.head(1), lambda q, t: (.5, 5))
            pd.DataFrame(dict(pair_id=['a'], query=['AC1'], text=['K'])).to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, 'ESM vocabulary'):
                plm_score.load_jobs(path)
        self.assertEqual(plm_score.average_precision([1, 0, 1], [.9, .8, .1]), (1 + 2/3)/2)
        self.assertEqual(plm_score.auroc([1, 0, 1], [.9, .8, .1]), .5)

    def test_score_designs(self):
        types = np.array(['diffusion']*4 + ['contact']*4)
        s = np.array([.9, .6, .8, .3, .5, .5, .7, .2])
        cases = research.score_designs(types, s, dict(permutation_mode='exhaustive', constant_control=True))
        self.assertEqual(sum(k.startswith('perm_diffusion') for k in cases), 23)
        self.assertEqual(sum(k.startswith('perm_contact') for k in cases), 11)   # tied scores
        np.testing.assert_allclose(cases['ppi_constant'][0], [.65]*4 + [.475]*4)
        for name, (values, kind) in cases.items():
            if kind:
                other = types != kind
                np.testing.assert_array_equal(values[other], s[other])  # other type untouched
                self.assertEqual(sorted(values[~other]), sorted(s[~other]))
        legacy = research.score_designs(types, s, dict(n_permutations=3))
        self.assertTrue(all(k in ['ppi'] or k.startswith('shuffled_') for k in legacy))

    def test_rank_metrics(self):
        rng = np.random.default_rng(0)
        receptor = rng.normal(size=200)
        target = receptor + rng.normal(scale=.3, size=200)
        signal = np.exp(receptor) + rng.normal(scale=.2, size=200)   # receptor + unrelated noise
        self.assertGreater(ev.spearman(signal, target), .7)
        self.assertLess(abs(ev.partial_spearman(signal, target, [receptor])), .15)
        # Signal fully determined by receptor expression carries no extra information.
        self.assertTrue(np.isnan(ev.partial_spearman(np.exp(receptor), target, [receptor])))
        self.assertEqual(ev.spearman(signal, target), ev.spearman(3*signal, target))  # scaling invariance
        self.assertTrue(np.isnan(ev.spearman(np.ones(5), np.arange(5))))

    def test_decomposition(self):
        base = np.ones((4, 4))
        blocks = [(slice(0, 2), slice(0, 2)), (slice(2, 4), slice(2, 4))]
        base[0:2, 2:4] = base[2:4, 0:2] = 0
        scaled = base.copy(); scaled[0:2, 0:2] *= 3
        result = analyze.decompose(base, scaled, blocks)
        self.assertGreater(result['between_pair_TV'], 0)
        self.assertAlmostEqual(result['within_pair_TV'], 0)
        moved = base.copy(); moved[0, 0], moved[1, 1] = 2, 0
        result = analyze.decompose(base, moved, blocks)
        self.assertAlmostEqual(result['between_pair_TV'], 0)
        self.assertGreater(result['within_pair_TV'], 0)

    def test_commot_mixture_and_pure_collective(self):
        with tempfile.TemporaryDirectory() as d, quiet():
            root = Path(d)
            research.demo(root/'in')
            commot_baseline.run(root/'in', root/'mix', dis_thr_h=(2,))
            commot_baseline.run(root/'in', root/'pure', dis_thr_h=(2,), weights=(1, 0, 0, 0))
            pure = pd.read_csv(root/'pure'/'summary.csv').iloc[0]
            self.assertLessEqual(pure.max_ligand_budget_excess, 1e-6)
            self.assertLessEqual(pure.max_receptor_budget_excess, 1e-6)
            with np.load(root/'mix'/'case_000.npz') as z:
                self.assertTrue((z['T'][~np.isfinite(z['C'])] == 0).all())
                self.assertGreater(np.trace(z['T'][:8, :8]), 0)  # same-spot transport allowed as in COMMOT
            self.assertTrue(json.loads((root/'pure'/'manifest.json').read_text())['budget_feasible'])


if __name__ == '__main__':
    unittest.main()
