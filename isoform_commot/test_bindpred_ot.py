"""Tests for bindpred_ot.py (reports/bindpred_ot_01/plan.md). Run: python -m unittest test_bindpred_ot
Needs vendor/_unot.py (OT solver), like test_stage2_s."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

import bindpred_ot as bo
import stage2_ppi as sp
import stage2_s as ss
from test_stage2_s import toy, run_args


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def proteins(tmp, inp):
    """Protein map (lr_symbol -> accession, sha256) and frozen sequences for the toy entities."""
    rng = np.random.default_rng(3)
    ents = pd.read_csv(inp/'entities.csv').entity.tolist()
    seqs = {e: ''.join(rng.choice(list('ACDEFGHIKLMNPQRSTVWY'), size=int(rng.integers(30, 60)))) for e in ents}
    pm = pd.DataFrame(dict(lr_symbol=ents, accession=['Q' + e for e in ents], length=[len(seqs[e]) for e in ents],
                           sequence_sha256=[sha(seqs[e]) for e in ents]))
    pm.to_csv(tmp/'pmap_full.csv', index=False)
    pd.DataFrame(dict(accession=pm.accession, taxon_id=10090, length=pm.length, sequence_sha256=pm.sequence_sha256,
                      uniprot_release='test', sequence=[seqs[e] for e in ents])).to_csv(tmp/'seq.tsv', sep='\t', index=False)
    return pm


def fake_predictions(jobs_dir, out, seed=0):
    pairs = pd.read_csv(jobs_dir/'pairs.csv')
    rng = np.random.default_rng(seed)
    pk = rng.uniform(6, 10, len(pairs))
    pd.DataFrame(dict(pair_id=pairs.pair_id, log10_kd_pred=-pk, pkd_hat=pk)).to_csv(out, index=False)


class Pipeline(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.inp = toy(self.tmp)
        self.pm = proteins(self.tmp, self.inp)
        bo.jobs(argparse.Namespace(input=[str(self.inp)], protein_map=str(self.tmp/'pmap_full.csv'),
                                   sequences=str(self.tmp/'seq.tsv'), output_dir=str(self.tmp/'jobs')))

    def test_jobs_orientation(self):
        ent = pd.read_csv(self.tmp/'jobs'/'entity_pairs.csv')
        pairs = pd.read_csv(self.tmp/'jobs'/'pairs.csv')
        self.assertEqual(len(ent), 5)
        for r in ent.itertuples():
            self.assertEqual(r.oriented_keys, f'Q{r.ligand}__Q{r.receptor}')     # ligand subunit first
        seq = pd.read_csv(self.tmp/'seq.tsv', sep='\t').set_index('accession').sequence
        for r in pairs.itertuples():
            self.assertEqual((r.lig_seq, r.rec_seq), (seq[r.lig_acc], seq[r.rec_acc]))

    def test_tampered_sequence_rejected(self):
        t = pd.read_csv(self.tmp/'seq.tsv', sep='\t')
        t.loc[0, 'sequence'] = t.loc[0, 'sequence'][::-1]
        t.to_csv(self.tmp/'bad.tsv', sep='\t', index=False)
        with self.assertRaises(ValueError):
            bo.jobs(argparse.Namespace(input=[str(self.inp)], protein_map=str(self.tmp/'pmap_full.csv'),
                                       sequences=str(self.tmp/'bad.tsv'), output_dir=str(self.tmp/'jobs2')))

    def test_import_best_pair_and_receptor_only(self):
        fake_predictions(self.tmp/'jobs', self.tmp/'bp.csv')
        bo.import_scores(argparse.Namespace(jobs_dir=str(self.tmp/'jobs'), bindpred=str(self.tmp/'bp.csv'),
                                            output_dir=str(self.tmp/'imp')))
        ppi = pd.read_csv(self.tmp/'imp'/'ppi_bindpred.csv')
        bp = pd.read_csv(self.tmp/'bp.csv').set_index('pair_id').pkd_hat
        for r in ppi.itertuples():
            self.assertAlmostEqual(r.pkd_hat, bp[f'Q{r.ligand}__Q{r.receptor}'])
            self.assertAlmostEqual(r.score, 10 ** (r.pkd_hat - bo.SCALE_PKD))
            self.assertEqual(r.best_subunit_pair, '__'.join(sorted([f'Q{r.ligand}', f'Q{r.receptor}'])))
        rec = pd.read_csv(self.tmp/'imp'/'ppi_bindpred_receptor_only.csv')
        for receptor, g in rec.groupby('receptor'):
            self.assertTrue(np.allclose(g.pkd_hat, ppi[ppi.receptor == receptor].pkd_hat.mean()))

    def test_end_to_end(self):
        fake_predictions(self.tmp/'jobs', self.tmp/'bp.csv')
        bo.import_scores(argparse.Namespace(jobs_dir=str(self.tmp/'jobs'), bindpred=str(self.tmp/'bp.csv'),
                                            output_dir=str(self.tmp/'imp')))
        for name in ('ppi_bindpred', 'ppi_bindpred_receptor_only'):
            sp.designs(argparse.Namespace(input=str(self.inp), ppi=str(self.tmp/'imp'/f'{name}.csv'),
                                          output=str(self.tmp/'imp'/f'designs_{name}_post.csv')))
        runs = self.tmp/'runs'
        bo.cases_cmd(argparse.Namespace(runs=str(runs), smoke=True))
        pmap = str(self.tmp/'pmap_full.csv')
        src = {'bp': ('designs_ppi_bindpred_post.csv', 'ppi_bindpred.csv'),
               'bprec': ('designs_ppi_bindpred_receptor_only_post.csv', 'ppi_bindpred_receptor_only.csv')}
        for sl in ('post', 'ant'):                       # the toy has one input; both slices use it
            for tag in ('bp', 'bprec', 'v12'):
                folder = runs/f'{sl}_{tag}'
                names = [c['name'] for c in json.loads((folder/'cases.json').read_text())['cases']]
                a = run_args(self.tmp, self.inp, folder/'cases.json', names, folder)
                if tag != 'v12':
                    a.designs, a.ppi_scores = (str(self.tmp/'imp'/f) for f in src[tag])
                    a.protein_map = pmap
                a.slice = sl
                ss.run(a)
        # experiment-2 stand-in: same cases written under the experiment-2 column names
        e2 = self.tmp/'exp2'; e2.mkdir()
        for sl in ('post', 'ant'):
            w, _ = ss.load_runs(runs/f'{sl}_bp'); v, _ = ss.load_runs(runs/f'{sl}_v12')
            pd.DataFrame(dict(no_ppi=w.no_ppi, snorm_b1=v.v12_snorm)).to_csv(e2/f'{sl}_main_receptor_partial_spearman.csv')
        bo.check_cmd(argparse.Namespace(runs=str(runs), exp2_summary=str(e2), output=str(self.tmp/'gate.json')))
        self.assertTrue(json.loads((self.tmp/'gate.json').read_text())['passed'])
        bad = pd.read_csv(e2/'post_main_receptor_partial_spearman.csv', index_col=0); bad['no_ppi'] += 1e-6
        bad.to_csv(e2/'post_main_receptor_partial_spearman.csv')
        with self.assertRaises(ValueError):
            bo.check_cmd(argparse.Namespace(runs=str(runs), exp2_summary=str(e2), output=str(self.tmp/'gate2.json')))
        bo.summarise(argparse.Namespace(runs=str(runs), output=str(self.tmp/'summary')))
        dec = json.loads((self.tmp/'summary'/'decisions.json').read_text())
        post = dec['details']['post']
        self.assertEqual(post['F1']['n_perm'], 2)
        self.assertIn(post['F2_vs_uniform']['label'], {'improved', 'worse', 'uncertain'})
        self.assertIn('bp_snorm_vs_bp_receptor_snorm', post['sensitivity'])
        self.assertIn(dec['verdict']['conclusion_post'], {'看不出 BindPred 的配对信息', '带有信息,但不超过平均分配',
                                                          'BindPred 的 S′ 改进了推断,且超过平均分配'})
        # the uniform and no_ppi cases do not depend on the scores: identical to the V12-designs runs
        u = ss.load_runs(runs/'post_bp')[0]
        self.assertFalse(np.allclose(u.bp_snorm.fillna(0), u.uniform.fillna(0)))


class Decisions(unittest.TestCase):
    def test_conclusion_table(self):
        self.assertEqual(bo.conclusion('information', 'improved'), 'BindPred 的 S′ 改进了推断,且超过平均分配')
        self.assertEqual(bo.conclusion('information', 'uncertain'), '带有信息,但不超过平均分配')
        self.assertEqual(bo.conclusion('no_information', 'improved'), '看不出 BindPred 的配对信息')
        self.assertEqual(bo.conclusion('undetermined', 'improved'), '无法判定(核心条件求解失败)')


if __name__ == '__main__':
    unittest.main()
