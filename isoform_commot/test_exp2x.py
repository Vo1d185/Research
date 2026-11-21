"""exp2x.py / exp2x_data.py (reports/stage2_ppi_03): length control, edges, pooled decisions, paper metric,
and equality with stage2_s.py on the toy input of test_stage2_s.py."""
import argparse
import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

import numpy as np
import pandas as pd

import exp2x as ex
import exp2x_data as xd
import stage2_s as ss
import test_stage2_s as ts


class LengthControl(unittest.TestCase):
    def pairs(self):
        return pd.DataFrame(dict(ligand=['A', 'A', 'A', 'B', 'C', 'C'], receptor=['X', 'Y', 'Z', 'X', 'Y', 'Z'],
                                 type='diffusion', s_raw=[0.9, 0.1, 0.5, 0.3, np.nan, 0.2],
                                 tokens=1.0, homophilic=False))

    def test_reassign_by_length(self):
        p = self.pairs()
        lens = dict(X=100.0, Y=900.0, Z=500.0)
        s, panel = ex.effective(p, dict(s='length', perm=0, override=None, panel='all'), lens)
        self.assertTrue(panel.all())
        np.testing.assert_allclose(s[:3], [0.1, 0.9, 0.5])          # longest Y gets 0.9, shortest X gets 0.1
        self.assertAlmostEqual(s[3], 0.3)                            # ligand with one receptor unchanged
        raw, _ = ss.effective(p, dict(s='raw', perm=0, override=None, panel='all'))
        self.assertAlmostEqual(s[4], raw[4])                         # unscored pair keeps the type median
        self.assertAlmostEqual(s[5], 0.2)                            # only one scored pair for C: unchanged
        for lig in 'ABC':                                            # per-ligand multiset preserved
            m = (p.ligand == lig).to_numpy()
            np.testing.assert_allclose(np.sort(s[m]), np.sort(raw[m]))

    def test_ties_by_name(self):
        p = self.pairs()
        s, _ = ex.effective(p, dict(s='length', perm=0, override=None, panel='all'), dict(X=5.0, Y=5.0, Z=5.0))
        np.testing.assert_allclose(s[:3], [0.9, 0.5, 0.1])           # equal lengths: X, Y, Z get descending values


class Edges(unittest.TestCase):
    def spatial(self):
        xy = np.array([[0, 0], [10, 0], [25, 0], [100, 0]], float)
        return pd.DataFrame(dict(x=xy[:, 0], y=xy[:, 1]), index=list('abcd'))

    def test_distance_rule(self):
        e = ex.edges_for(self.spatial(), dict(contact_rule='distance', contact_cutoff_um=12.0), 30.0)
        s, t, d = e['contact']
        got = sorted(zip(s.tolist(), t.tolist()))
        self.assertEqual(got, [(0, 0), (0, 1), (1, 0), (1, 1), (2, 2), (3, 3)])
        s, t, d = e['diffusion']
        self.assertEqual(sorted(zip(s.tolist(), t.tolist())),
                         [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (3, 3)])
        self.assertTrue((d <= 30.0).all())

    def test_none_rule_and_hex(self):
        import stage2
        e = ex.edges_for(self.spatial(), dict(contact_rule='none'), 30.0)
        self.assertEqual(e['contact'][0].tolist(), [0, 1, 2, 3])
        tmp = Path(tempfile.mkdtemp()); inp = ts.toy(tmp)
        sp = pd.read_csv(inp/'spatial.csv', index_col=0)
        a, b = ex.edges_for(sp, {}, 363.0757), dict(zip(['diffusion', 'contact'], stage2.spatial_edges(sp, 363.0757)))
        for k in a:
            for u, v in zip(a[k], b[k]):
                np.testing.assert_array_equal(u, v)


class Pooled(unittest.TestCase):
    def test_paired_labels(self):
        c1 = pd.DataFrame(dict(raw=[0.2, 0.3, 0.25], no_ppi=[0.1, 0.1, 0.1]))
        c2 = pd.DataFrame(dict(raw=[0.2, 0.25], no_ppi=[0.1, 0.1]))
        r = ex.pooled_paired({'breast': c1, 'seqfish': c2}, 'raw', 'no_ppi')
        self.assertAlmostEqual(r['delta'], (0.1 + 0.2 + 0.15 + 0.1 + 0.15) / 5)
        self.assertEqual(r['label'], 'improved')
        c2b = pd.DataFrame(dict(raw=[0.0, 0.05], no_ppi=[0.1, 0.1]))   # one dataset negative: not 'improved'
        self.assertNotEqual(ex.pooled_paired({'breast': c1, 'seqfish': c2b}, 'raw', 'no_ppi')['label'], 'improved')

    def test_pooled_d1(self):
        rng = np.random.default_rng(0)
        res = {}
        for ds, shift in [('breast', 0.05), ('seqfish', 0.05)]:
            cols = {'raw': [0.1 + shift] * 4}
            for k in range(1, 21):
                cols[f'raw_perm_{k:03d}'] = list(0.1 + 0.001 * rng.standard_normal(4))
            common = pd.DataFrame(cols)
            res[ds] = (common, {}, common)
        r = ex.pooled_d1(res, 'raw')
        self.assertEqual(r['n_perm'], 20)
        self.assertAlmostEqual(r['p'], 1 / 21)
        self.assertEqual(r['label'], 'information')


class PaperMetric(unittest.TestCase):
    def test_nan_rules(self):
        E = np.array([[1, 0], [2, 0], [3, 0], [4, 0]], float)
        items = [('L', 'R1', ['g1']), ('L', 'R2', ['g9']), ('L', 'R3', ['g1'])]
        sig = {('L', 'R1'): np.array([1, 2, 3, 4.0]), ('L', 'R2'): np.array([1, 2, 3, 4.0])}
        m = ex.paper_metric(sig, E, ['g1', 'g2'], items)
        self.assertAlmostEqual(m.spearman[0], 1.0)
        self.assertEqual(m.spearman[1], 0.0)          # no target among the evaluation genes -> 0
        self.assertEqual(m.spearman[2], 0.0)          # pair without signal (absent) -> zeros -> 0
        self.assertTrue(m.signal_zero[2])

    def test_paper_items(self):
        tmp = Path(tempfile.mkdtemp())
        base = 'downstream_gene_correlation/scseqcomm/data/'
        with zipfile.ZipFile(tmp/'z.zip', 'w') as z:
            z.writestr(base + 'TF_PPR_REACTOME_mouse.csv', pd.DataFrame(dict(
                receptor=['Ra,Rb', 'Rc', 'Rc', 'Rd'], tf=['T1', 'T2', 'T3', 'T4'],
                tf_PPR=[1e-3, 1e-3, 1e-6, 1e-3])).to_csv())
            z.writestr(base + 'TF_list-TF_TG_TRRUSTv2_RegNetwork_High_mouse.txt', pd.DataFrame(dict(x=['T1', 'T2', 'T3'])).to_csv(index=False))
            for tf, tg in [('T1', ['g1', 'g2']), ('T2', ['g3']), ('T3', ['g4'])]:
                z.writestr(f'{base}TF_TG_TRRUSTv2_RegNetwork_High_mouse/{tf}.csv', pd.DataFrame(dict(x=tg)).to_csv(index=False))
        lr = pd.DataFrame(dict(ligand=['L1', 'L2', 'L3', 'L4'], receptor=['Ra_Rb', 'Rb_Ra', 'Rc', 'Rd']))
        items = ex.paper_items(lr, 'mouse', tmp/'z.zip')
        self.assertEqual([(l, r) for l, r, _ in items], [('L1', 'Ra_Rb'), ('L3', 'Rc')])   # order-sensitive match;
        self.assertEqual(items[0][2], ['g1', 'g2'])                                       # Rd: TF without list
        self.assertEqual(items[1][2], ['g3'])                                             # T3 below 1e-5


class Cases(unittest.TestCase):
    def test_counts(self):
        full, smoke = ex.case_list(), ex.case_list(smoke=True)
        self.assertEqual(len(full), 5 + 2 * ex.N_PERM)
        self.assertEqual(len(smoke), 9)
        self.assertTrue(all(c['beta'] == (0.0 if c['s'] == 'none' else 1.0) for c in full))
        self.assertEqual(len({c['name'] for c in full}), len(full))


class SameAsStage2S(unittest.TestCase):
    """exp2x.run must give the same receptor metrics as stage2_s.run for no_ppi and raw (beta 1) on a hex input."""

    def test_equal_metrics(self):
        tmp = Path(tempfile.mkdtemp()); inp = ts.toy(tmp)
        ss_cases = tmp/'ss_cases.json'
        ss.cases_cmd(argparse.Namespace(setting='main', slice='post', output=str(ss_cases), smoke=True))
        ss.run(ts.run_args(tmp, inp, ss_cases, ['no_ppi', 'raw'], tmp/'ss_out'))
        ex_cases = tmp/'ex_cases.json'
        ex.cases_cmd(argparse.Namespace(dataset='toy', output=str(ex_cases), smoke=True))
        a = ts.run_args(tmp, inp, ex_cases, ['no_ppi', 'raw', 'length_raw'], tmp/'ex_out')
        pm = pd.read_csv(tmp/'pmap.csv').assign(lr_symbol=['A', 'X', 'Y', 'C', 'Z'])   # real maps carry lr_symbol
        pm.to_csv(tmp/'pmap_sym.csv', index=False)
        a.protein_map, a.dataset = str(tmp/'pmap_sym.csv'), 'toy'
        ex.run(a)
        for name in ['no_ppi', 'raw']:
            x = pd.read_csv(tmp/'ss_out'/name/'metrics.csv'); y = pd.read_csv(tmp/'ex_out'/name/'metrics.csv')
            pd.testing.assert_frame_equal(x, y)
        self.assertTrue((tmp/'ex_out'/'length_raw'/'metrics.csv').exists())
        ref = pd.DataFrame({n: pd.read_csv(tmp/'ss_out'/n/'metrics.csv').query('top_tf == 10').set_index('receptor')
                            .partial_spearman for n in ['no_ppi', 'raw']})
        ref.to_csv(tmp/'ref.csv')
        ex.regression(argparse.Namespace(runs=str(tmp/'ex_out'), reference=str(tmp/'ref.csv'),
                                         output=str(tmp/'reg.json'), case=['no_ppi', 'raw']))
        self.assertTrue(json.loads((tmp/'reg.json').read_text())['passed'])


class DataConstants(unittest.TestCase):
    def test_paper_table(self):
        self.assertEqual(set(xd.PAPER), {'post', 'breast', 'seqfish'})
        self.assertEqual(xd.PAPER['seqfish']['input_scale'], 'log1p')
        self.assertAlmostEqual(xd.PAPER['post']['median'], 0.180)
        self.assertEqual(xd.SEQFISH_UM_PER_PX, 0.103)


if __name__ == '__main__':
    unittest.main()
