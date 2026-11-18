"""gpu_unot.py must reproduce the official COMMOT solver. Runs on CPU (torch device 'cpu'); on the server the
same tests also run with device 'cuda' when available. Skipped if torch or the vendored solver is missing."""
from pathlib import Path
import unittest

import numpy as np

try:
    import torch
except ImportError:          # commot-research has no torch: these tests run in the GPU environment
    torch = None

HAVE_SOLVER = (Path(__file__).parent/'vendor'/'_unot.py').exists()
DEVICES = ['cpu'] + (['cuda'] if torch is not None and torch.cuda.is_available() else [])


def problem(seed=0, ns=60, nt=50, density=0.3):
    rng = np.random.default_rng(seed)
    mask = rng.random((ns, nt)) < density
    mask[np.arange(ns), rng.integers(0, nt, ns)] = True
    mask[rng.integers(0, ns, nt), np.arange(nt)] = True
    r, s = np.nonzero(mask)
    c = rng.uniform(0, 3, len(r))
    a, b = rng.random(ns), rng.random(nt)
    return a, b, r, s, c


@unittest.skipUnless(torch is not None and HAVE_SOLVER, 'needs torch and vendor/_unot.py')
class Port(unittest.TestCase):
    def test_fixed_iterations_match(self):
        import gpu_unot as gu
        for dev in DEVICES:
            for eps, rho in ((0.1, 10.0), (0.2, 1.0)):
                a, b, r, s, c = problem()
                res = gu.compare(*gu.normalised_problem(r, s, c, a, b), eps, rho, 300, dev)
                self.assertTrue(res['order_preserved'])
                self.assertLess(res['max_rel_diff'], 1e-12 if dev == 'cpu' else gu.EQUIVALENCE_TOL)
                self.assertEqual(res['iterations'], 301)

    def test_full_solve_matches_research(self):
        import gpu_unot as gu
        from research import solve_sparse
        for dev in DEVICES:
            for eps, rho in ((0.1, 10.0), (0.2, 1.0)):
                a, b, r, s, c = problem(seed=1)
                T0, i0 = solve_sparse(r, s, c, a, b, eps, rho, 200000, 1e-11)
                T1, i1 = gu.solve_sparse_gpu(r, s, c, a, b, eps, rho, 200000, 1e-11, device=dev)
                self.assertEqual(i0['niter'], i1['niter'])
                self.assertEqual(i0['cap_reached'], i1['cap_reached'])
                self.assertLess(np.abs(T0 - T1).max() / np.abs(T0).max(), 1e-9)

    def test_cap_is_reported(self):
        import gpu_unot as gu
        a, b, r, s, c = gu.normalised_problem(*problem(seed=2)[2:], *problem(seed=2)[:2])
        _, n = gu.unot_sinkhorn_l1_sparse_torch(a, b, r, s, c, 0.1, 10.0, nitermax=20, stopthr=0.0, device='cpu')
        self.assertEqual(n, 21)


if __name__ == '__main__':
    unittest.main()
