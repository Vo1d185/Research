import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from cellchat_inventory_grid import Library, entity_detection, gene_fractions, genes_by_key, load_visium, run_grid


class GridTests(unittest.TestCase):
    def test_complex_detection_codetect_vs_per_subunit(self):
        # 3 spots; A in spots 0,1; B in spots 1,2 -> co-detected only in spot 1
        X = sparse.csc_matrix(np.array([[1, 0], [2, 3], [0, 1]]))
        genes = pd.Index(['A', 'B'])
        det = entity_detection({'A_B': ['A', 'B'], 'A': ['A'], 'A_C': ['A', 'C']},
                               X, genes, gene_fractions(X, 3))
        self.assertAlmostEqual(det.loc['A_B', 'codetect'], 1/3)
        self.assertAlmostEqual(det.loc['A_B', 'per_subunit'], 2/3)
        self.assertAlmostEqual(det.loc['A', 'codetect'], 2/3)
        self.assertTrue(np.isnan(det.loc['A_C', 'codetect']))  # missing subunit
        self.assertEqual(det.loc['A_C', 'missing'], 'C')

    def test_visium_duplicate_symbols_are_summed(self):
        import h5py
        M = sparse.csc_matrix(np.array([[1, 0], [0, 2], [3, 0]]))  # features x barcodes
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'m.h5'
            with h5py.File(path, 'w') as f:
                g = f.create_group('matrix')
                g['data'], g['indices'], g['indptr'], g['shape'] = M.data, M.indices, M.indptr, M.shape
                feats = g.create_group('features')
                feats['name'] = np.array([b'Ccl21a', b'Ccl21a', b'Gfap'])
                feats['feature_type'] = np.array([b'Gene Expression'] * 3)
            v = load_visium(path)
        self.assertEqual(list(v['genes']), ['Ccl21a', 'Gfap'])
        np.testing.assert_array_equal(v['X'].toarray(), [[1, 3], [2, 0]])
        self.assertEqual(v['duplicated_symbols'], ['Ccl21a'])

    def test_library_any_all_upper_and_ranking(self):
        rt = pd.DataFrame(dict(receptor=['FZD4', 'FZD4', 'FZD4'], tf=['T1', 'T2', 'T3'], weight=[.1, .5, .5]))
        tt = pd.DataFrame(dict(tf=['T2', 'T3'], target=['GENEX', 'GENEY']))
        lib = Library('ccc', rt, tt, upper=True)
        linked, tfs = lib.receptor(['Fzd4', 'Lrp6'])
        self.assertEqual(linked, ['Fzd4'])            # any: linked; all: not (Lrp6 missing)
        self.assertEqual(tfs, ['T2', 'T3', 'T1'])     # weight desc, ties by name
        keymap = genes_by_key(pd.Index(['Genex', 'Geney', 'Fzd4']), upper=True)
        self.assertEqual(lib.target_symbols(['T2', 'T3'], keymap), {'Genex', 'Geney'})

    def test_single_gene_rule_matches_evaluate_targets(self):
        from evaluate_targets import target_sets
        rng = np.random.default_rng(1)
        genes = [f'G{i}' for i in range(40)] + ['L', 'R']
        X = sparse.csc_matrix(rng.poisson(.3, (50, len(genes))))
        frac = gene_fractions(X, 50)
        rt = pd.DataFrame(dict(receptor='R', tf=[f'T{i}' for i in range(12)], weight=rng.uniform(size=12)))
        tt = pd.DataFrame([(f'T{i}', f'G{(3*i + k) % 40}') for i in range(12) for k in range(4)], columns=['tf', 'target'])
        ref = target_sets(['R'], {'R': 'R'}, rt, tt, dict(zip(genes, frac)), {'L', 'R'}, 10, 5, .05).iloc[0]
        pairs = pd.DataFrame(dict(interaction_name=['L_R'], annotation=['Secreted Signaling'], pathway_name=['P'],
                                  ligand=['L'], receptor=['R'], ligand_complex=[False], receptor_complex=[False]))
        X[:, -2:] = 1  # L and R detected everywhere
        ds = {'toy': dict(X=sparse.csc_matrix(X), genes=pd.Index(genes), n_spots=50, duplicated_symbols=[])}
        steps, passing = run_grid(pairs, {'L': ['L'], 'R': ['R']}, ds, [Library('sc', rt, tt, upper=False)],
                                  thresholds=(.05,))
        got = passing[passing.complex_rule == 'single_gene_only'].iloc[0]
        self.assertEqual(got.n_targets, ref.n_targets)
        self.assertEqual(bool(got.step_c), ref.status == 'evaluated')


if __name__ == '__main__':
    unittest.main()
