import io
from pathlib import Path
import unittest
import zipfile

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from simulation import benchmark as bm
from simulation import fields
from simulation.pde import Grid, Model, System

ARCHIVE = Path('data/zenodo_7272562/pde_model_comparison.zip')


def blob(centre, width=12.0):
    xy = fields.grid_points()
    f = np.exp(-0.5 * (((xy[..., 0] - centre[0]) ** 2 + (xy[..., 1] - centre[1]) ** 2) / width ** 2))
    f[f < 0.1] = 0
    return f


def small_scenario():
    """Two ligands, one receptor (pairs L1-R1, L2-R1), compact supports: a small OT problem on the 50 x 50 grid."""
    p = np.stack([blob((-20, 0)), blob((20, 0))])
    r = np.stack([blob((0, 0), 16.0)])
    return dict(p=p, r0=r, binding=np.array([[1.0], [1.0]]), structure=np.array(0), repeat=np.array(1),
                source=np.array('test'))


class Fields(unittest.TestCase):
    def test_reproducible_and_shape_of_fields(self):
        params = fields.FieldParams(5, 40.0, 0.25)
        a = fields.scenario_fields(3, 1, 2, 3, params)
        b = fields.scenario_fields(3, 1, 2, 3, params)
        c = fields.scenario_fields(3, 2, 2, 3, params)
        for x, y in zip(a, b):
            np.testing.assert_array_equal(x, y)
        self.assertFalse(np.array_equal(a[0], c[0]))
        for f in np.concatenate(a):
            self.assertLessEqual(f.max(), 1.0)
            self.assertGreaterEqual(f[f > 0].min(), fields.THRESHOLD)
            self.assertGreater(f.max(), 0.5)                      # a blob centre lies within the domain

    def test_statistics(self):
        f = np.zeros((10, 10)); f[1:3, 1:3] = 0.5; f[6:9, 6:9] = 1.0
        self.assertAlmostEqual(fields.support(f), 13 / 100)
        self.assertEqual(fields.blobs(f), 2)
        diag = np.zeros((10, 10)); diag[0, 0] = diag[1, 1] = 1        # diagonal neighbours: 4-connectivity -> 2
        self.assertEqual(fields.blobs(diag), 2)
        g = np.zeros((10, 10)); g[1:3, 1:3] = 1.0                       # receptor support of 4 points
        self.assertAlmostEqual(fields.overlap(f, g), 1.0)
        self.assertAlmostEqual(fields.overlap(np.roll(f, 1, axis=0), g), 0.5)

    def test_calibration_picks_minimum_with_grid_order_ties(self):
        params = fields.FieldParams(3, 40.0, 0.5)
        binding = np.ones((2, 2))
        reference = fields.pooled([fields.scenario_statistics(*fields.scenario_fields(9, s, 2, 2, params), binding)
                                   for s in range(1, 40)])
        old_grid, old_rep = fields.CALIBRATION_GRID, fields.CALIBRATION_REPEATS
        fields.CALIBRATION_GRID = dict(k_max=[1, 3, 8], sigma_hi=[20.0, 40.0], sigma_ratio=[0.5])
        fields.CALIBRATION_REPEATS = range(1001, 1021)
        try:
            rows, chosen = fields.calibrate(reference, {9: binding})
        finally:
            fields.CALIBRATION_GRID, fields.CALIBRATION_REPEATS = old_grid, old_rep
        objective = [r['objective'] for r in rows]
        self.assertEqual(len(rows), 6)
        best = rows[objective.index(min(objective))]
        self.assertEqual((best['k_max'], best['sigma_hi'], best['sigma_ratio']),
                         (chosen.k_max, chosen.sigma_hi, chosen.sigma_ratio))
        self.assertEqual((chosen.k_max, chosen.sigma_hi), (3, 40.0))


class Truth(unittest.TestCase):
    def test_jacobian_matches_finite_differences(self):
        grid = Grid(6, -50, 50, 2)
        rng = np.random.default_rng(3)
        model = Model(D=[200, 50], a=[[1, 1, 0], [0, 1, 1]], b=[[1e-3, 1e-3, 1e-3], [1e-3, 1e-3, 1e-3]], c=[0.2, 0.1])
        s = System(model, grid, rng.uniform(size=(2, 6, 6)), rng.uniform(size=(3, 6, 6)))
        y = rng.uniform(size=s.nfields * s.size)
        J = s.jacobian(0, y).toarray()
        h, fd = 1e-6, np.empty_like(J)
        for c in range(len(y)):
            e = np.zeros_like(y); e[c] = h
            fd[:, c] = (s.rhs(0, y + e) - s.rhs(0, y - e)) / (2 * h)
        np.testing.assert_allclose(J, fd, atol=1e-6)

    def test_steady_checks_and_dynamic_convergence(self):
        scen = small_scenario()
        sub = dict(scen, p=scen['p'][:, 18:32, 18:32], r0=scen['r0'][:, 18:32, 18:32])
        model = Model(D=[50.0, 50.0], a=[[1.0], [1.0]], b=[[1e-3], [1e-3]], c=[0.2, 0.2])
        system = System(model, Grid(14, -56, 56, 2), sub['p'], sub['r0'])
        L, R, C = system.steady_state()
        info = bm.steady_checks(system, L, R, C)
        self.assertTrue(info['passed'], info)
        (Ld, Rd, Cd), t, history = system.run_to_steady_state(t_chunk=1000.0, tol=1e-8, t_max=1e6, method='BDF',
                                                              rtol=1e-10, atol=1e-14, per_field=True)
        for q in range(len(C)):
            self.assertLessEqual(np.abs(Cd[q] - C[q]).max() / C[q].max(), bm.CHECKS['dynamic_match'])

    def test_system_for_uses_lambda_grid(self):
        scen = small_scenario()
        for k in bm.KS:
            s = bm.system_for(scen, k)
            self.assertAlmostEqual(np.sqrt(s.m.D[0] / s.m.c[0]), bm.lam(k))
            self.assertEqual(s.grid.dx, bm.H)
        self.assertAlmostEqual(bm.lam(0) / bm.H, 3.953, places=3)


class Conditions(unittest.TestCase):
    def test_names_and_sharing(self):
        self.assertEqual(bm.run_name(0, 'cd_2lambda', 0.1, 10, 50000), bm.run_name(1, 'cd_1lambda', 0.1, 10, 50000))
        self.assertEqual(bm.run_name(-1, 'cd_4lambda', 0.1, 10, 50000), bm.run_name(1, 'cd_1lambda', 0.1, 10, 50000))
        self.assertNotEqual(bm.run_name(0, 'r0_2', 0.1, 10, 50000), bm.run_name(1, 'r0_1', 0.1, 10, 50000))
        self.assertEqual(bm.parse_condition('r0_0.125', 2), ('r0', 0.125 * bm.lam(2)))
        self.assertEqual(bm.parse_condition('cd_150', -2), ('cd', 150.0))
        per_k = pd.Series([k for k, _ in bm.job_conditions()]).value_counts()
        self.assertEqual(per_k[-2], 11)                    # cd at 1 lambda (< h) skipped
        self.assertTrue((per_k.drop(-2) == 12).all())

    def test_entries_match_stage2_construction(self):
        import stage2
        scen = small_scenario()
        cutoff, r0 = 4 * bm.lam(-1), 0.5 * bm.lam(-1)
        row, col, d, pair = bm.ot_entries(scen['p'], scen['r0'], scen['binding'], cutoff)
        xy = fields.grid_points().reshape(-1, 2)
        n = len(xy)
        pairs = cKDTree(xy).query_pairs(cutoff, output_type='ndarray')
        s = np.concatenate([np.arange(n), pairs[:, 0], pairs[:, 1]]); t = np.concatenate([np.arange(n), pairs[:, 1], pairs[:, 0]])
        dist = np.linalg.norm(xy[s] - xy[t], axis=1)
        lr = pd.DataFrame(dict(ligand=['L1', 'L2'], receptor=['R1', 'R1'], type=['diffusion', 'diffusion']))
        a, b = scen['p'].reshape(2, -1).reshape(-1), scen['r0'].reshape(1, -1).reshape(-1)
        r2, c2, cost2, _, _ = stage2.build_sparse_cost(lr, ['L1', 'L2'], ['R1'], a, b, n,
                                                       {'diffusion': (s, t, dist)}, r0)
        mine = sorted(zip(row.tolist(), col.tolist(), np.round(np.log1p(d / r0), 12).tolist()))
        theirs = sorted(zip(r2.tolist(), c2.tolist(), np.round(cost2, 12).tolist()))
        self.assertEqual(mine, theirs)

    def test_ours_and_commot_runs(self):
        scen = small_scenario()
        sig, info = bm.solve_ours(scen, 'r0', bm.lam(0), 4 * bm.lam(0), 0.2, 1.0, 50000, 1e-11)
        self.assertTrue(info['check_passed'], info)
        self.assertEqual(info['forbidden_mass'], 0.0)
        self.assertEqual(sig.shape, (2, 2500))
        self.assertLessEqual(sig.sum(axis=0).max(), scen['r0'].max() + 1e-9)      # receptor budget per point
        self.assertAlmostEqual(sig.sum(), info['total_T'])
        sig, info = bm.solve_ours(scen, 'inv', None, 4 * bm.lam(0), 0.2, 1.0, 50000, 1e-11)
        self.assertTrue(info['check_passed'], info)
        sig_c, info_c = bm.solve_commot(scen, 4 * bm.lam(0), 0.2, 1.0, 50000, 1e-11)
        self.assertTrue(info_c['check_passed'], info_c)
        self.assertEqual(info_c['forbidden_mass'], 0.0)
        self.assertAlmostEqual(sig_c.sum(), info_c['total_T'])
        sig_r, info_r = bm.solve_commot(scen, 150.0, 0.25, 1.0, 10000, 1e-8, square=True)
        self.assertEqual(sig_r.shape, (2, 2500))


class Metrics(unittest.TestCase):
    def test_pair_metrics_rules(self):
        r0 = np.array([[0.0, 1.0], [1.0, 1.0]])
        truth = np.array([[9.0, 1.0], [2.0, 3.0]])
        m = bm.pair_metrics(truth, np.array([[100.0, 0.5], [1.0, 2.0]]), r0)        # point (0,0) excluded
        self.assertAlmostEqual(m['spearman'], 1.0)
        self.assertAlmostEqual(m['rmse'], np.sqrt(np.mean([0.25, 1.0, 1.0])))
        m = bm.pair_metrics(truth, np.ones((2, 2)), r0)                              # constant signal: NaN -> 0
        self.assertEqual((m['spearman'], m['spearman_nan']), (0.0, True))

    @unittest.skipUnless(ARCHIVE.exists(), 'archive not downloaded')
    def test_metrics_reproduce_archive_evaluation_rule(self):
        """Archive evaluate_result_*.py on case 2 sim 1 (archived PDE result at the last time, archived COT output)."""
        import h5py
        from sklearn.metrics import mean_squared_error
        zf = zipfile.ZipFile(ARCHIVE)
        base = 'pde_model_comparison/data/case_2/sim_1/'
        cot = bm.h5ad_matrix(zf.read(base + 'cot/adata.h5ad'))
        with h5py.File(io.BytesIO(zf.read(base + 'pde/result.h5')), 'r') as h:
            data = h['data'][:]
        labels = zf.read(base + 'pde/result_labels.txt').decode().split()
        for i, j in [(0, 0), (0, 1), (1, 1), (1, 2)]:
            rec = np.load(io.BytesIO(zf.read(base + f'pde/receptor_init_{j+1}.npy')))
            pde = data[-1, labels.index(f'Complex_{i+1}_{j+1}')]
            sig = cot[f'receiver-L{i+1}-R{j+1}'].to_numpy(float)
            idx = np.where(rec.flatten() > 0)[0]
            r, _ = spearmanr(pde.flatten()[idx], sig[idx])
            r = 0 if np.isnan(r) else r
            rmse = np.sqrt(mean_squared_error(pde.flatten()[idx], sig[idx]))
            mine = bm.pair_metrics(pde, sig.reshape(pde.shape), rec)
            self.assertEqual(mine['spearman'], r)
            self.assertEqual(mine['rmse'], rmse)


if __name__ == '__main__':
    unittest.main()
