import unittest

import numpy as np

from simulation.pde import Grid, Model, System, bump, crescent, gaussian2d


def two_ligands_one_receptor(n=80, D=(200.0, 200.0), l1=300.0, l2=100.0):
    grid = Grid(n, 0.0, 400.0, 1)
    x = grid.coords()
    model = Model(D=D, a=[[1.0], [1.0]], b=[[1e-3], [1e-3]], c=[0.2, 0.2])
    return System(model, grid, [bump(x, l1, 50), bump(x, l2, 50)], [bump(x, 200, 50)])


class Operators(unittest.TestCase):
    def test_laplacian_is_zero_flux_and_matches_matrix(self):
        rng = np.random.default_rng(1)
        for grid in (Grid(30, -1, 1, 1), Grid(12, -1, 1, 2)):
            u = rng.uniform(size=grid.shape)
            lap = grid.laplacian(u)
            self.assertAlmostEqual(lap.sum(), 0.0, places=8)                      # no flux through the boundary
            np.testing.assert_allclose(grid.laplacian_matrix() @ u.ravel(), lap.ravel(), rtol=1e-12, atol=1e-9)

    def test_receptor_plus_complex_conserved_pointwise(self):
        grid = Grid(10, -50, 50, 2)
        rng = np.random.default_rng(2)
        model = Model(D=[200, 50], a=[[1, 1, 0], [0, 1, 1]], b=[[1e-3, 1e-3, 0], [0, 1e-3, 1e-3]], c=[0.2, 0.1])
        s = System(model, grid, rng.uniform(size=(2, 10, 10)), rng.uniform(size=(3, 10, 10)))
        y = rng.uniform(size=s.nfields*s.size)
        _, dR, dC = s.split(s.rhs(0, y))
        for j in range(3):
            total = dR[j] + sum(dC[k] for k, (i, jj) in enumerate(model.pairs) if jj == j)
            np.testing.assert_allclose(total, 0, atol=1e-12)


class SteadyState(unittest.TestCase):
    def test_closed_form_is_a_fixed_point(self):
        s = two_ligands_one_receptor()
        L, R, C = s.steady_state()
        residual = s.rhs(0, np.concatenate([L.ravel(), R.ravel(), C.ravel()]))
        self.assertLess(np.abs(residual).max(), 1e-9)
        np.testing.assert_allclose(R + C.sum(axis=0), s.r0, rtol=1e-12)          # receptor conserved

    def test_long_integration_reaches_closed_form(self):
        s = two_ligands_one_receptor(n=40)
        (L, R, C), t, _ = s.run_to_steady_state(t_chunk=2000.0, tol=1e-9)
        L0, R0, C0 = s.steady_state()
        np.testing.assert_allclose(C, C0, rtol=1e-5, atol=1e-8 * C0.max())
        np.testing.assert_allclose(L, L0, rtol=1e-5, atol=1e-8 * L0.max())

    def test_symmetric_case_is_mirrored_and_biased_case_favours_the_closer_ligand(self):
        L, R, C = two_ligands_one_receptor().steady_state()
        np.testing.assert_allclose(C[0], C[1][::-1], rtol=1e-8, atol=1e-12)
        _, _, C = two_ligands_one_receptor(l1=270.0).steady_state()
        self.assertGreater(C[0].sum(), C[1].sum())
        _, _, C = two_ligands_one_receptor(D=(200.0, 800.0)).steady_state()
        self.assertGreater(C[1].sum(), C[0].sum())                               # faster-diffusing ligand 2


class Shapes(unittest.TestCase):
    def test_shapes(self):
        g = Grid(50, -200, 200, 2).coords()
        self.assertAlmostEqual(gaussian2d(g, (4, 4), (20, 10), 30).max(), 1.0, places=2)
        arc = crescent(g, (0, 0), 150, 15, -100, 100)
        self.assertGreater(arc[g[..., 0] > 100].sum(), 10 * arc[g[..., 0] < -100].sum())   # right side only
        self.assertEqual(bump(np.array([0.0, 50.0, 100.0]), 0.0, 50.0).tolist(), [1.0, 0.0, 0.0])


if __name__ == '__main__':
    unittest.main()
