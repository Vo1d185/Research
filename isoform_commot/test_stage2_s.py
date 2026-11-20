"""stage2_s.py (reports/stage2_ppi_02): score handling, evaluation equal to evaluate_targets.py, decisions."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

import stage2_s as ss

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def toy(tmp, cutoff=363.0757):
    """Stage-2 style input (16 spots, 5 pairs), gene export, scSeqComm-like db, designs, ppi scores, protein map."""
    tmp = Path(tmp)
    inp = tmp/'in'; inp.mkdir()
    pts = [(r, c) for r in range(4) for c in range(8) if (r + c) % 2 == 0]
    spots = [f's{i}' for i in range(len(pts))]
    pd.DataFrame(dict(x=[c * 50.0 for r, c in pts], y=[r * 86.6025 for r, c in pts], array_row=[r for r, _ in pts],
                      array_col=[c for _, c in pts]), index=spots).to_csv(inp/'spatial.csv', index_label='spot')
    lr = pd.DataFrame(dict(ligand=['A', 'A', 'B', 'C', 'C'], receptor=['X', 'Y', 'X', 'Y', 'Z'],
                           type=['diffusion', 'diffusion', 'diffusion', 'contact', 'contact']))
    lr.to_csv(inp/'lr.csv', index=False)
    ents = ['A', 'B', 'C', 'X', 'Y', 'Z']
    rng = np.random.default_rng(0)
    expr = pd.DataFrame({e: np.round(rng.random(len(spots)) * (rng.random(len(spots)) > 0.2), 4) for e in ents}, index=spots)
    expr.to_csv(inp/'expression.csv', index_label='spot')
    pd.DataFrame(dict(entity=ents, subunits=ents)).to_csv(inp/'entities.csv', index=False)
    (inp/'config.json').write_text(json.dumps(dict(eps=0.2, rho=1.0, diffusion_cutoff_um=cutoff, nitermax=5000, stopthr=1e-9)))
    files = ['expression.csv', 'spatial.csv', 'lr.csv', 'entities.csv', 'config.json']
    (inp/'preparation.json').write_text(json.dumps(dict(status='complete', neighbour_spacing_um=100.0,
                                                        prepared_sha256={f: sha(inp/f) for f in files})))
    genes = ents + [f'G{i}' for i in range(1, 13)]
    X = np.c_[expr.to_numpy(), rng.random((len(spots), 12)) + 0.01]
    g = tmp/'genes'; g.mkdir()
    np.savez_compressed(g/'gene_expression.npz', X=X, spots=np.asarray(spots, str), genes=np.asarray(genes, str))
    (g/'gene_export.json').write_text(json.dumps(dict(status='complete', gene_expression_sha256=sha(g/'gene_expression.npz'))))
    pd.DataFrame(dict(lr_symbol=ents + [None] * 12, data_symbol=genes, detected_fraction=1.0)).to_csv(g/'genes.csv', index=False)
    db = tmp/'db'; db.mkdir()
    pd.DataFrame(dict(receptor=['X', 'Y', 'Z', 'X'], tf=['T1', 'T2', 'T3', 'T2'], weight=[0.9, 0.8, 0.7, 0.1],
                      source='KEGG')).to_csv(db/'receptor_tf.csv', index=False)
    pd.DataFrame(dict(tf=['T1'] * 6 + ['T2'] * 6 + ['T3'] * 6,
                      target=[f'G{i}' for i in list(range(1, 7)) + list(range(4, 10)) + list(range(7, 13))])).to_csv(db/'tf_tg.csv', index=False)
    (db/'db.json').write_text(json.dumps(dict(status='complete', converted_sha256={f: sha(db/f) for f in ['receptor_tf.csv', 'tf_tg.csv']})))
    designs = lr.assign(s_raw=[0.9, 0.05, np.nan, 0.6, 0.001])
    designs.to_csv(tmp/'designs.csv', index=False)
    (tmp/'designs.csv.json').write_text(json.dumps(dict(input_sha256={'lr.csv': sha(inp/'lr.csv')})))
    pd.DataFrame(dict(ligand=lr.ligand, receptor=lr.receptor, type=lr.type, homophilic=[False, False, False, True, False],
                      best_subunit_pair=['PA__PX', 'PA__PY', None, 'PC__PY', 'PC__PZ'])).to_csv(tmp/'ppi.csv', index=False)
    pd.DataFrame(dict(accession=['PA', 'PX', 'PY', 'PC', 'PZ'], length=[300, 1500, 200, 100, 400])).to_csv(tmp/'pmap.csv', index=False)
    return inp


def run_args(tmp, inp, cases, names, out, save=None, setting='main'):
    return argparse.Namespace(input=str(inp), designs=str(tmp/'designs.csv'), ppi_scores=str(tmp/'ppi.csv'),
                              protein_map=str(tmp/'pmap.csv'), genes=str(tmp/'genes'), db=str(tmp/'db'), cases=str(cases),
                              output=str(out), setting=setting, slice='post', case=names, device='cpu', save_case=save)


class Scores(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp()); toy(self.tmp)
        lr = pd.read_csv(self.tmp/'in'/'lr.csv')
        self.pairs = ss.pair_table(lr, pd.read_csv(self.tmp/'designs.csv'), pd.read_csv(self.tmp/'ppi.csv'),
                                   pd.read_csv(self.tmp/'pmap.csv'))

    def spec(self, **kw):
        return {**dict(s='raw', perm=0, override=None, panel='all'), **kw}

    def test_impute_normalise_uniform(self):
        s, panel = ss.effective(self.pairs, self.spec())
        self.assertTrue(panel.all())
        self.assertAlmostEqual(s[2], np.median([0.9, 0.05]))                  # unscored diffusion pair: type median
        s, _ = ss.effective(self.pairs, self.spec(s='snorm'))
        self.assertAlmostEqual(s[0] + s[1], 1.0); self.assertAlmostEqual(s[2], 1.0); self.assertAlmostEqual(s[3] + s[4], 1.0)
        s, _ = ss.effective(self.pairs, self.spec(s='uniform'))
        self.assertTrue(np.allclose(s, [0.5, 0.5, 1.0, 0.5, 0.5]))
        s, panel = ss.effective(self.pairs, self.spec(panel='scored'))
        self.assertEqual(int(panel.sum()), 4); self.assertEqual(len(s), 4)

    def test_permutation_and_overrides(self):
        base, _ = ss.effective(self.pairs, self.spec())
        for k in (1, 2, 3):
            s, _ = ss.effective(self.pairs, self.spec(perm=k))
            self.assertEqual(sorted(s[[0, 1]]), sorted(base[[0, 1]]))        # diffusion scores only among themselves
            self.assertEqual(sorted(s[[3, 4]]), sorted(base[[3, 4]]))
            self.assertAlmostEqual(s[2], base[2])                              # imputed value never permuted
        s, _ = ss.effective(self.pairs, self.spec(override='long'))           # PA__PX = 1803 tokens > 1603
        self.assertAlmostEqual(s[0], np.median([0.9, 0.05])); self.assertAlmostEqual(s[1], 0.05)
        s, _ = ss.effective(self.pairs, self.spec(override='homo'))
        self.assertAlmostEqual(s[3], np.median([0.6, 0.001])); self.assertAlmostEqual(s[4], 0.001)

    def test_case_lists(self):
        main = ss.case_list('main', 'post')
        self.assertEqual(sum(c['group'] == 'perm' for c in main), 200)
        self.assertEqual({c['name'] for c in main if c['group'] == 'core'}, set(ss.CORE))
        self.assertEqual([c['beta'] for c in main if c['name'] in ('raw', 'snorm', 'uniform')], [1.0, 0.2, 0.2])
        second = ss.case_list('second', 'post')
        self.assertEqual({c['name'] for c in second}, set(ss.CORE) | {'raw_n200k', 'snorm_n200k'})
        self.assertEqual(len(ss.case_list('second', 'ant')), 4)


class EndToEnd(unittest.TestCase):
    def test_run_equivalence_summarise(self):
        import evaluate_targets as et
        tmp = Path(tempfile.mkdtemp()); inp = toy(tmp)
        runs = tmp/'runs'/'post_main'
        ss.cases_cmd(argparse.Namespace(setting='main', slice='post', output=str(runs/'cases.json'), smoke=True))
        names = [c['name'] for c in json.loads((runs/'cases.json').read_text())['cases']]
        ss.run(run_args(tmp, inp, runs/'cases.json', ['no_ppi', 'raw'], runs, save=str(tmp/'full')))
        ss.run(run_args(tmp, inp, runs/'cases.json', names, runs))             # resumable: no_ppi, raw skipped
        self.assertTrue(all((runs/n/'case.json').exists() for n in names))
        self.assertTrue((runs/'raw'/'signal.npz').exists()); self.assertFalse((runs/'raw_perm_001'/'signal.npz').exists())
        et.evaluate(argparse.Namespace(input=str(inp), genes=str(tmp/'genes'), db=str(tmp/'db'),
                                       results=[str(tmp/'full'/'no_ppi'), str(tmp/'full'/'raw')], output=str(tmp/'eval'),
                                       top_tf=[10, 5, 20], target_score=['mean', 'zmean'], min_targets=5, min_detected=0.05,
                                       covariate=None, commot_rule=False, complex_link='any'))
        ss.equivalence(argparse.Namespace(runs=str(runs), eval_dir=str(tmp/'eval'), case=['no_ppi', 'raw'],
                                          output=str(tmp/'equivalence.json')))
        self.assertTrue(json.loads((tmp/'equivalence.json').read_text())['passed'])
        raw = pd.read_csv(runs/'raw'/'metrics.csv'); base = pd.read_csv(runs/'no_ppi'/'metrics.csv')
        self.assertEqual(set(raw.receptor), {'X', 'Y', 'Z'})
        self.assertFalse(np.allclose(raw.total_received, base.total_received))  # S changed the allocation
        ss.summarise(argparse.Namespace(runs=str(tmp/'runs'), output=str(tmp/'summary')))
        dec = json.loads((tmp/'summary'/'decisions.json').read_text())
        self.assertIn('raw', dec['main']['verdict']); self.assertEqual(dec['main']['details']['post']['D1_raw']['n_perm'], 2)


class Cutoff(unittest.TestCase):
    """Amendment 01: the plan's cutoff is used even when config.json stores another one (stage2_input_post stores 500)."""
    def test_plan_cutoff_used(self):
        import stage2
        import stage2_inject as si
        res = {}
        for cut in (363.0757, 500.0):
            tmp = Path(tempfile.mkdtemp()); inp = toy(tmp, cut)
            _, cfg, _, spatial, *_, edges = ss.load_input(inp)
            self.assertEqual(cfg['diffusion_cutoff_um'], ss.CUTOFF_UM); self.assertEqual(cfg['config_cutoff_um'], cut)
            ref = stage2.spatial_edges(spatial, ss.CUTOFF_UM)[0]
            for u, v in zip(edges['diffusion'], ref):
                np.testing.assert_array_equal(u, v)
            if cut == 500.0:   # the toy is large enough for 500 to add edges, so the test can fail
                self.assertGreater(len(si.load(inp)[-1]['diffusion'][0]), len(ref[0]))
            runs = tmp/'runs'/'post_main'
            ss.cases_cmd(argparse.Namespace(setting='main', slice='post', output=str(runs/'cases.json'), smoke=True))
            ss.run(run_args(tmp, inp, runs/'cases.json', ['no_ppi', 'raw'], runs))
            rec = json.loads((runs/'raw'/'case.json').read_text())
            self.assertEqual((rec['cutoff_um'], rec['config_cutoff_um']), (ss.CUTOFF_UM, cut))
            res[cut] = {n: pd.read_csv(runs/n/'metrics.csv') for n in ['no_ppi', 'raw']}
        for n in ['no_ppi', 'raw']:
            pd.testing.assert_frame_equal(res[363.0757][n], res[500.0][n])


class Decisions(unittest.TestCase):
    def test_paired_and_d1(self):
        rng = np.random.default_rng(0)
        base = pd.Series(rng.normal(size=40), index=[f'r{i}' for i in range(40)])
        wide = pd.DataFrame(dict(no_ppi=base, raw=base + 0.05 + rng.normal(scale=0.01, size=40), snorm=base - 0.05,
                                 uniform=base + rng.normal(scale=0.01, size=40)))
        for k in range(1, 20):
            wide[f'raw_perm_{k:03d}'] = base + rng.normal(scale=0.01, size=40)
            wide[f'snorm_perm_{k:03d}'] = base + rng.normal(scale=0.01, size=40)
        info = {n: dict(group='perm' if 'perm' in n else 'core', error='', panel='all') for n in wide}
        d = ss.slice_decisions(wide, info, 'main')
        self.assertEqual(d['D1_raw']['label'], 'information'); self.assertAlmostEqual(d['D1_raw']['p'], 1 / 20)
        self.assertEqual(d['D1_snorm']['label'], 'no_information')
        self.assertEqual(d['D2_raw']['label'], 'improved'); self.assertEqual(d['D2_snorm']['label'], 'worse')
        self.assertEqual(d['D3_raw_vs_snorm']['label'], 'raw_better')
        self.assertEqual(ss.conclusion('information', 'uncertain'), '带有配对信息,但目前不足以改进推断')
        self.assertEqual(ss.conclusion('no_information', 'improved'), '在本评价下看不出配对信息')

    def test_stability(self):
        tmp = Path(tempfile.mkdtemp())
        sig = np.abs(np.random.default_rng(1).normal(size=(3, 50)))
        for name, s in [('raw', sig), ('raw_n200k', sig * (1 + 1e-6)), ('snorm', sig), ('snorm_n200k', sig * 1.01)]:
            (tmp/name).mkdir()
            np.savez_compressed(tmp/name/'signal.npz', receptors=np.array(['a', 'b', 'c']), signal=s)
            (tmp/name/'case.json').write_text(json.dumps(dict(total_transport=float(s.sum()))))
        r = ss.stability(tmp)
        self.assertTrue(r['raw']['passed']); self.assertFalse(r['snorm']['passed'])


if __name__ == '__main__':
    unittest.main()
