"""Toy tests for stage2_inject.py: truth formulas, score calibration, normalisation, decision rules."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

import stage2_inject as si


class Truth(unittest.TestCase):
    def test_competitive_binding(self):
        x = pd.DataFrame(dict(L1=[1.0, 0.0], L2=[0.0, 1.0], R=[2.0, 2.0]))
        lr = pd.DataFrame(dict(ligand=['L1', 'L2'], receptor=['R', 'R'], type=['contact', 'contact']))
        L = {'L1': np.array([1.0, 0.5]), 'L2': np.array([0.0, 2.0])}
        C, load = si.complexes(lr, np.array([1.0, 3.0]), L, x)
        # spot 1: load = 1*0.5 + 3*2 = 6.5; C1 = 2*0.5/7.5, C2 = 2*6/7.5
        np.testing.assert_allclose(C[:, 1], [1.0 / 7.5, 12.0 / 7.5])
        np.testing.assert_allclose(C[:, 0], [2 * 1 / 2, 0.0])
        self.assertAlmostEqual(load['R'][1], 6.5)

    def test_ligand_fields(self):
        x = pd.DataFrame(dict(A=[1.0, 0.0, 0.0], B=[0.0, 2.0, 0.0]))
        lr = pd.DataFrame(dict(ligand=['A', 'B'], receptor=['R', 'R'], type=['diffusion', 'contact']))
        dif = (np.array([0, 1, 0]), np.array([0, 1, 1]), np.array([0.0, 0.0, 100.0]))
        con = (np.array([1, 1]), np.array([1, 2]), np.array([0.0, 100.0]))
        L = si.ligand_fields(x, lr, {'diffusion': dif, 'contact': con}, lam=100.0)
        np.testing.assert_allclose(L['A'], [1.0, np.exp(-1.0), 0.0])
        np.testing.assert_allclose(L['B'], [0.0, 2.0, 2.0])


class Scores(unittest.TestCase):
    def test_calibration_and_latent(self):
        rng = np.random.default_rng(0)
        K = 10 ** rng.uniform(-1, 1, 200); decoy = rng.random(200) < 1 / 3
        u = si.latent(K, decoy)
        self.assertTrue((u[decoy] == -1).all() and (u[~decoy] >= 0).all() and (u[~decoy] <= 2).all())
        for target in (0.9, 0.7, 0.6):
            sig = si.sigma_for(u, decoy, target)
            self.assertAlmostEqual(si.expected_auroc(u, decoy, sig), target, places=4)
        self.assertEqual(si.sigma_for(u, decoy, 1.0), 0.0)

    def test_ligand_normalisation(self):
        lr = pd.DataFrame(dict(ligand=['A', 'A', 'B'], receptor=['X', 'Y', 'X']))
        s = si.effective_score(lr, [0.6, 0.2, 0.3], 'ligand')
        np.testing.assert_allclose(s, [0.75, 0.25, 1.0])
        np.testing.assert_allclose(si.effective_score(lr, [0.0, 0.5, 1.0], 'raw'), [si.S_FLOOR, 0.5, 1.0])

    def test_received(self):
        T = np.array([1.0, 2.0, 3.0]); col = np.array([0 * 2 + 1, 1 * 2 + 1, 1 * 2 + 0]); p = np.array([0, 1, 1])
        np.testing.assert_allclose(si.received(T, col, p, 2, 2), [[0, 1], [3, 2]])


def table(rows):
    return pd.DataFrame(rows, columns=['seed', 'variant', 'norm', 'beta', 'pooled_spearman',
                                       'receptor_spearman_mean', 'cap_reached', 'eps', 'decoy_transport_fraction'])


class Decide(unittest.TestCase):
    def grid(self, oracle_gain, perm_gain=0.0, other_gain=0.0, cap=False):
        rows = []
        for seed in (1, 2, 3):
            rows.append((seed, 'none', 'none', 0.0, 0.5, 0.5, False, 0.1, 0.3))
            for b in (0.1, 1.0):
                rows.append((seed, 'uniform', 'uniform', b, 0.5, 0.5, False, 0.1, 0.3))
                for v in si.AUROC_LEVELS:
                    g = oracle_gain if v == 'q100' else other_gain
                    for norm in ('raw', 'ligand'):
                        rows.append((seed, v, norm, b, 0.5 + g * b, 0.5 + g * b, cap and b == 1.0, 0.1, 0.2))
                        rows.append((seed, v + '_perm', norm, b, 0.5 + perm_gain, 0.5, False, 0.1, 0.3))
        return table(rows)

    def test_positive(self):
        dec, _ = si.decide(self.grid(0.2), dict(auroc=0.85, ci95=[0.8, 0.9]))
        self.assertTrue(dec['R1_oracle_usable']); self.assertTrue(dec['R5_receptor_metric_sees_oracle'])
        self.assertEqual(dec['R3_threshold_auroc'], 1.0)
        self.assertEqual(dec['R4_v12']['placement'], 'below')

    def test_null_and_convergence(self):
        dec, _ = si.decide(self.grid(0.0))
        self.assertFalse(dec['R1_oracle_usable']); self.assertIsNone(dec['R3_threshold_auroc'])
        dec, _ = si.decide(self.grid(0.2, cap=True))
        self.assertTrue(dec['R1_cells']['raw|0.1']); self.assertFalse(dec['R1_cells']['raw|1'])
        self.assertEqual(dec['runs_not_converged'], 3 * 5 * 2)

    def test_threshold_levels(self):
        dec, _ = si.decide(self.grid(0.2, other_gain=0.1))
        self.assertEqual(dec['R3_threshold_auroc'], 0.6); self.assertTrue(dec['R3_monotone'])


@unittest.skipUnless((Path(__file__).parent/'vendor'/'_unot.py').exists(), 'run fetch_solver.py first')
class EndToEnd(unittest.TestCase):
    def test_tiny(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); inp = tmp/'in'; inp.mkdir()
            pts = [(r, c) for r in range(4) for c in range(8) if (r + c) % 2 == 0]
            spots = [f's{i}' for i in range(len(pts))]
            pd.DataFrame(dict(x=[c * 50.0 for r, c in pts], y=[r * 86.6025 for r, c in pts],
                              array_row=[r for r, _ in pts], array_col=[c for _, c in pts]), index=spots
                         ).to_csv(inp/'spatial.csv', index_label='spot')
            lr = pd.DataFrame(dict(ligand=['A', 'A', 'B', 'C', 'C'], receptor=['X', 'Y', 'X', 'Y', 'Z'],
                                   type=['diffusion', 'diffusion', 'diffusion', 'contact', 'contact']))
            lr.to_csv(inp/'lr.csv', index=False)
            ents = ['A', 'B', 'C', 'X', 'Y', 'Z']
            rng = np.random.default_rng(0)
            pd.DataFrame({e: np.round(rng.random(len(spots)) * (rng.random(len(spots)) > 0.3), 4) for e in ents},
                         index=spots).to_csv(inp/'expression.csv', index_label='spot')
            pd.DataFrame(dict(entity=ents, subunits=ents)).to_csv(inp/'entities.csv', index=False)
            (inp/'config.json').write_text(json.dumps(dict(eps=0.2, rho=1.0, diffusion_cutoff_um=363.0757,
                                                          nitermax=5000, stopthr=1e-9)))
            files = ['expression.csv', 'spatial.csv', 'lr.csv', 'entities.csv', 'config.json']
            (inp/'preparation.json').write_text(json.dumps(dict(status='complete', prepared_sha256={
                f: hashlib.sha256((inp/f).read_bytes()).hexdigest() for f in files})))
            si.truth(argparse.Namespace(input=inp, output=tmp/'t', seed=1, lambda_um=100.0, decoy_fraction=0.4,
                                        log10k_range=2.0))
            si.scores(argparse.Namespace(truth=tmp/'t', output=tmp/'s', seed=1))
            ns = dict(input=inp, truth=tmp/'t', scores=tmp/'s', seed=1, beta=[0.2], r0_um=50.0, eps=None, rho=None,
                      diffusion_cutoff_um=None, stopthr=None, nitermax=None)
            si.run(argparse.Namespace(output=tmp/'r1', variant='none', norm='raw', **ns))
            si.run(argparse.Namespace(output=tmp/'r2', variant='q100', norm='ligand', **ns))
            a = pd.read_csv(tmp/'r1'/'metrics.csv'); b = pd.read_csv(tmp/'r2'/'metrics.csv')
            self.assertEqual(len(a), 1); self.assertEqual(len(b), 1)
            self.assertTrue(np.isfinite(b.pooled_spearman).all())
            self.assertLessEqual(b.decoy_transport_fraction.iloc[0], a.decoy_transport_fraction.iloc[0] + 1e-9)


if __name__ == '__main__':
    unittest.main()
