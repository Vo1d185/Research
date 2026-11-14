import argparse
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

import plans
import stage2
from research import backend_sparse, solve, solve_sparse, sparse_kkt_check, sparse_mass_check


def hex_slide(rows=6, cols=12):
    """Visium-like lattice: same-row neighbours 2 columns apart; 100 um spacing."""
    pts = [(r, c) for r in range(rows) for c in range(cols) if (r + c) % 2 == 0]
    ids = [f's{r}_{c}' for r, c in pts]
    return pd.DataFrame(dict(x=[c*50.0 for r, c in pts], y=[r*50.0*np.sqrt(3) for r, c in pts],
                             array_row=[r for r, c in pts], array_col=[c for r, c in pts]), index=ids)


def write_input(root, spatial, expression, lr, cfg):
    root.mkdir(parents=True)
    expression.to_csv(root/'expression.csv', index_label='spot')
    spatial.to_csv(root/'spatial.csv', index_label='spot')
    lr.to_csv(root/'lr.csv', index=False)
    ents = sorted(set(lr.ligand) | set(lr.receptor))
    pd.DataFrame(dict(entity=ents, subunits=ents)).to_csv(root/'entities.csv', index=False)
    (root/'config.json').write_text(json.dumps(cfg))
    files = ['expression.csv', 'spatial.csv', 'lr.csv', 'entities.csv', 'config.json']
    (root/'preparation.json').write_text(json.dumps(dict(status='complete', sample='synthetic',
        neighbour_spacing_um=100.0, prepared_sha256={f: stage2.digest(root/f) for f in files})))


class SparseSolver(unittest.TestCase):
    def test_matches_dense_elementwise(self):
        rng = np.random.default_rng(3)
        for eps, rho in [(0.2, 1.0), (0.1, 10.0), (0.5, 0.5)]:
            C = rng.uniform(0, 2, (30, 25)); C[rng.uniform(size=C.shape) < .6] = np.inf
            C[rng.uniform(size=C.shape) < .05] = 0.0       # explicit zero costs must stay allowed
            a, b = rng.uniform(0, 1, 30), rng.uniform(0, 1, 25); a[:2] = 0; b[-1] = 0
            Td = solve(C, a, b, eps, rho)
            r, c = np.nonzero(np.isfinite(C))
            Ts, info = solve_sparse(r, c, C[r, c], a, b, eps, rho)
            D = np.zeros_like(Td); D[r, c] = Ts
            self.assertLessEqual(np.abs(D - Td).max(), 1e-10 * Td.max())
            self.assertEqual(Td[~np.isfinite(C)].sum(), 0)
            self.assertEqual(info['kkt_rows'], info['active_rows'])     # every row and column covered
            self.assertLessEqual(info['kkt_tol'], 1e-9)

    def test_mass_check_detects_wrong_plans(self):
        rng = np.random.default_rng(4)
        C = rng.uniform(0, 1, (12, 10))
        a, b = rng.uniform(.1, 1, 12), rng.uniform(.1, 1, 10)
        r, c = np.nonzero(np.isfinite(C))
        T, info = solve_sparse(r, c, C[r, c], a, b, 0.2, 1.0)
        self.assertEqual(info['rows_checked'], 12)
        s = max(a.sum(), b.sum())
        bad = T / s
        bad[r == 3] *= 1.05                                 # one row off its optimum
        with self.assertRaises(ValueError):
            sparse_mass_check(r, c, C[r, c], bad, a/s, b/s, 0.2, 1.0)
        with self.assertRaises(ValueError):                 # transport computed with another cost
            sparse_mass_check(r, c, C[r, c] * 1.1, T / s, a/s, b/s, 0.2, 1.0)
        for args in [(r, c, C[r, c], bad, a/s, b/s, 0.2, 1.0), (r, c, C[r, c] * 1.1, T / s, a/s, b/s, 0.2, 1.0)]:
            with self.assertRaises(ValueError):
                sparse_kkt_check(*args)

    def test_kkt_check_covers_rows_the_primal_check_cannot(self):
        rng = np.random.default_rng(8)
        C = rng.uniform(0, 1, (40, 30)); C[rng.uniform(size=C.shape) < .7] = np.inf
        a, b = rng.uniform(.1, 1, 40), rng.uniform(.1, 1, 30)
        r, c = np.nonzero(np.isfinite(C))
        T, info = solve_sparse(r, c, C[r, c], a, b, 0.1, 10.0)    # COMMOT's eps, rho
        self.assertLess(info['rows_checked'], info['rows'])         # unused mass below double precision
        self.assertEqual(info['kkt_rows'], 40)
        s = max(a.sum(), b.sum())
        coo = sparse.coo_matrix((C[r, c], (r, c)), shape=C.shape)
        early = backend_sparse()(a/s, b/s, coo, 0.1, 10.0, nitermax=3, stopthr=1e-11)   # not converged
        with self.assertRaisesRegex(ValueError, 'not converged'):
            sparse_kkt_check(early.row, early.col, early.data * 0 + C[early.row, early.col], early.data,
                             a/s, b/s, 0.1, 10.0)

    def test_empty_and_duplicates(self):
        T, info = solve_sparse([], [], [], [1., 2.], [1.], 0.2, 1.0)
        self.assertEqual(len(T), 0)
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            solve_sparse([0, 0], [0, 0], [0., 1.], [1.], [1.], 0.2, 1.0)


class CommotPureCollective(unittest.TestCase):
    def test_cot_sparse_equals_combine_with_weights_1000(self):
        import commot_baseline
        combine, collective = commot_baseline.official_cot()
        rng = np.random.default_rng(7)
        xy = rng.uniform(0, 5, (15, 2))
        M = np.hypot(*(xy[:, None] - xy[None]).transpose(2, 0, 1))
        S, D = rng.uniform(0, 1, (15, 2)) * (rng.uniform(size=(15, 2)) > .3), rng.uniform(0, 1, (15, 3))
        A = np.array([[1.0, np.inf, 1.0], [1.0, 1.0, np.inf]])
        cut = 2.5 * np.ones_like(A)
        P1 = collective(S, D, A, M, cut, eps_p=0.1, rho=10.0, nitermax=50000, stopthr=1e-11)
        P2 = combine(S, D, A, M, cut, eps_p=0.1, rho=10.0, weights=(1.0, 0.0, 0.0, 0.0), nitermax=50000,
                     stopthr=1e-11)
        self.assertEqual(set(P1), set(P2))
        for key in P1:
            np.testing.assert_array_equal(P1[key].toarray(), P2[key].toarray())


class Geometry(unittest.TestCase):
    def test_hex_neighbours_and_edges(self):
        sp = hex_slide()
        s, t = stage2.hex_pairs(sp.array_row.to_numpy(), sp.array_col.to_numpy())
        d = np.hypot(*(sp[['x', 'y']].to_numpy()[s] - sp[['x', 'y']].to_numpy()[t]).T)
        np.testing.assert_allclose(d, 100.0)
        counts = np.bincount(s, minlength=len(sp))
        self.assertEqual(counts.max(), 6)                    # interior spots
        self.assertLess(counts.min(), 6)                     # edge spots have fewer neighbours
        (ds, dt, dd), (cs, ct, cd) = stage2.spatial_edges(sp, 250.0)
        self.assertTrue(((ds == dt) == (dd == 0)).all())     # self edges included, d = 0
        self.assertTrue((dd <= 250).all())
        self.assertEqual(len(cs), len(sp) + len(s))          # self + hex neighbours

    def test_sparse_cost_equals_independent_dense_cost(self):
        sp = hex_slide(4, 8)
        n = len(sp)
        lr = pd.DataFrame(dict(ligand=['L1', 'L1', 'L2'], receptor=['R1', 'R2', 'R1'],
                               type=['diffusion', 'contact', 'contact']))
        ligands, receptors = ['L1', 'L2'], ['R1', 'R2']
        rng = np.random.default_rng(5)
        a = rng.uniform(0, 1, 2*n) * (rng.uniform(size=2*n) > .2)
        b = rng.uniform(0, 1, 2*n) * (rng.uniform(size=2*n) > .2)
        edges = dict(zip(['diffusion', 'contact'], stage2.spatial_edges(sp, 300.0)))
        row, col, c, d, p = stage2.build_sparse_cost(lr, ligands, receptors, a, b, n, edges, 120.0)
        C = stage2.dense_cost(sp, lr, ligands, receptors, 300.0, 120.0)
        usable = np.isfinite(C) & (a[:, None] > 0) & (b[None, :] > 0)
        self.assertEqual(set(zip(*np.nonzero(usable))), set(zip(row.tolist(), col.tolist())))
        np.testing.assert_allclose(C[row, col], c, rtol=0, atol=1e-14)
        self.assertTrue((c[lr.type.to_numpy()[p] == 'contact'] == 0).all())


class EndToEnd(unittest.TestCase):
    def test_run_compare_dense_and_received_signal(self):
        sp = hex_slide(5, 10)
        rng = np.random.default_rng(6)
        x = pd.DataFrame(rng.poisson(2, (len(sp), 4)).astype(float), index=sp.index, columns=['L1', 'L2', 'R1', 'R2'])
        lr = pd.DataFrame(dict(ligand=['L1', 'L2', 'L1'], receptor=['R1', 'R1', 'R2'],
                               type=['diffusion', 'contact', 'contact']))
        cfg = dict(stage2.DEFAULTS, r0_um=[100.0, 200.0], diffusion_cutoff_um=300.0,
                   coordinate_unit='micrometer', expression_level='gene', expression_state='processed')
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            write_input(root/'in', sp, x, lr, cfg)
            # command-line values differ from config.json; compare-dense must use the values the run used
            stage2.run(argparse.Namespace(input=str(root/'in'), output=str(root/'out'), r0_um=None,
                                          eps=0.2, rho=1.0, diffusion_cutoff_um=250.0))
            stage2.compare_dense(argparse.Namespace(input=str(root/'in'), results=str(root/'out')))
            table = pd.read_csv(root/'out'/'dense_comparison.csv')
            self.assertTrue(table.passed.all() and table.same_eligible_pattern.all())
            self.assertLessEqual(table.max_rel_T_difference.max(), 1e-6)
            plan = plans.load(root/'out'/'case_000.npz')
            self.assertTrue(plan['sparse'])
            dense = dict(plan, sparse=False, T=plans.to_coo(plan).toarray())
            sp_pair, sp_rec = plans.received_signal(plan, lr, len(sp))
            de_pair, de_rec = plans.received_signal(dense, lr, len(sp))
            for k in sp_pair:
                np.testing.assert_allclose(sp_pair[k], de_pair[k], rtol=1e-12)
            for k in sp_rec:
                np.testing.assert_allclose(sp_rec[k], de_rec[k], rtol=1e-12)
            summary = pd.read_csv(root/'out'/'summary.csv')
            self.assertEqual((summary.eps.iloc[0], summary.rho.iloc[0], summary.diffusion_cutoff_um.iloc[0]),
                             (0.2, 1.0, 250.0))
            np.testing.assert_allclose(summary.frac_ligand, summary.total_T / x[['L1', 'L2']].to_numpy().sum())
            manifest = json.loads((root/'out'/'manifest.json').read_text())
            self.assertEqual(manifest['overridden'], dict(eps=0.1, rho=10.0, diffusion_cutoff_um=300.0))
            self.assertLessEqual(plan['D'].max(), 250.0)
            self.assertGreater(summary.total_T_contact.iloc[0], 0)
            self.assertEqual(summary.check_rows_checked.iloc[0], summary.check_rows.iloc[0])


class Report(unittest.TestCase):
    def test_r0_rule_tie_and_undefined_receptors(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            cut = stage2.COMMOT_TUTORIAL_CUTOFF_UM
            ours, commot, ev = root/'A', root/'B', root/'eval'
            for folder in (ours, commot, ev):
                folder.mkdir()
            (ours/'manifest.json').write_text(json.dumps(dict(method='research sparse (stage 2)',
                config=dict(nitermax=50000, stopthr=1e-11))))
            pd.DataFrame(dict(case=['case_000', 'case_001', 'case_002'], r0_um=[50.0, 100.0, 200.0], eps=0.1,
                              rho=10.0, diffusion_cutoff_um=cut, frac_ligand=.5, frac_receptor=.4, total_T=1.0,
                              mean_distance=10.0, solve_seconds=1.0)).to_csv(ours/'summary.csv', index=False)
            (commot/'manifest.json').write_text(json.dumps(dict(method='COMMOT original')))
            pd.DataFrame(dict(case=['case_000'], r0=[cut], dis_thr_um=[cut], eps=0.1, rho=10.0, scheme='collective',
                              nitermax=50000, stopthr=1e-11, frac_ligand=.5, frac_receptor=.4, total_T=1.0,
                              mean_distance=10.0, solve_seconds=1.0, check_passed=True)).to_csv(commot/'summary.csv',
                                                                                              index=False)
            pd.DataFrame(dict(receptor=['R1', 'R2', 'R3'], top_tf=10, status='evaluated')).to_csv(ev/'target_sets.csv',
                                                                                                  index=False)
            # R1, R2 defined at every r0: means 0.3, 0.3, 0.2 -> tie between 50 and 100, rule picks 50.
            # R3 undefined at r0 = 200: excluded from the selection although it favours 100.
            partial = {('A', 'case_000'): [.2, .4, .0], ('A', 'case_001'): [.4, .2, .9],
                       ('A', 'case_002'): [.2, .2, np.nan], ('B', 'case_000'): [.1, .1, .1]}
            rows = [dict(results=str(root/f), case=c, level='receptor', top_tf=10, target_score='mean',
                         score_mode='no_ppi', receptor=r, partial_spearman=v, spearman=v, receptor_only_spearman=0.0)
                    for (f, c), values in partial.items() for r, v in zip(['R1', 'R2', 'R3'], values)]
            pd.DataFrame(rows).to_csv(ev/'case_metrics.csv', index=False)
            (ev/'eval_manifest.json').write_text(json.dumps(dict(status='complete', results=[str(ours), str(commot)],
                                                                 primary_top_tf=10, target_score=['mean'])))
            stage2.report(argparse.Namespace(eval=str(ev), output=str(root/'report')))
            out = json.loads((root/'report'/'report.json').read_text())
            self.assertEqual(out['selected_r0_um'], 50.0)
            self.assertEqual(out['receptors_for_r0_selection'], 2)
            paired = pd.read_csv(root/'report'/'paired_receptors.csv').iloc[0]
            self.assertEqual((paired.ours, paired.baseline, paired.receptors), ('A', 'B', 3))
            self.assertAlmostEqual(paired.mean_difference, np.mean([.1, .3, -.1]))


class Preprocessing(unittest.TestCase):
    def test_normalize_filter_and_complex_minimum(self):
        X = sparse.csr_matrix(np.array([[1, 0, 3, 0], [0, 2, 2, 0], [4, 0, 0, 5]], float))
        Y, target = stage2.normalize_log1p(X)
        self.assertEqual(target, 4.0)                       # median of totals 4, 4, 9
        np.testing.assert_allclose(Y.toarray()[2], np.log1p(np.array([4, 0, 0, 5]) * 4/9))
        var = pd.Index(['A', 'B', 'C', 'D'])
        frac = np.asarray((X > 0).sum(axis=0)).ravel() / 3   # 2/3, 1/3, 2/3, 1/3
        pairs = pd.DataFrame(dict(ligand=['A', 'A', 'A'], receptor=['B_C', 'C_D', 'E'], type='diffusion',
                                  pathway_name='p'))
        subunits = {'A': ['A'], 'B_C': ['B', 'C'], 'C_D': ['C', 'D'], 'E': ['E']}
        audit = stage2.commot_filter(pairs, subunits, var, frac, 0.5)
        self.assertEqual(list(audit.reason), ['subunit_below_min_cell_pct', 'subunit_below_min_cell_pct',
                                              'subunit_not_measured'])
        audit = stage2.commot_filter(pairs, subunits, var, frac, 0.3)
        self.assertEqual(list(audit.reason[:2]), ['selected', 'selected'])  # each subunit separately
        e = stage2.entity_expression(Y, var, ['B_C'], subunits)
        np.testing.assert_allclose(e['B_C'], np.minimum(Y.toarray()[:, 1], Y.toarray()[:, 2]))

    def test_complex_receptor_target_sets(self):
        from evaluate_targets import target_sets
        rt = pd.DataFrame(dict(receptor=['Fzd4', 'Fzd4', 'Lrp6'], tf=['T1', 'T2', 'T2'], weight=[.2, .1, .9]))
        tt = pd.DataFrame(dict(tf=['T1']*3 + ['T2']*3, target=[f'G{i}' for i in range(6)]))
        detected = {f'G{i}': 1.0 for i in range(6)}
        symbol = {s: s for s in ['Fzd4', 'Lrp6', 'Ror1']}
        sets = target_sets(['Fzd4_Lrp6'], symbol, rt, tt, detected, set(), 1, 1, .05,
                           subunits={'Fzd4_Lrp6': ['Fzd4', 'Lrp6']})
        self.assertEqual(sets.iloc[0].tfs, 'T2')            # T2 = max(0.1, 0.9) over subunits
        sets = target_sets(['Fzd4_Ror1'], symbol, rt, tt, detected, set(), 2, 1, .05,
                           subunits={'Fzd4_Ror1': ['Fzd4', 'Ror1']})
        self.assertEqual(sets.iloc[0].status, 'evaluated')  # any subunit linked
        sets = target_sets(['Fzd4_Ror1'], symbol, rt, tt, detected, set(), 2, 1, .05,
                           subunits={'Fzd4_Ror1': ['Fzd4', 'Ror1']}, link_rule='all')
        self.assertEqual(sets.iloc[0].status, 'not_all_subunits_linked')


if __name__ == '__main__':
    unittest.main()
