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


class Amendment02(unittest.TestCase):
    """amendment_02.md: convergence = solver checks passed (Stage 2 section 8.7); stability between two caps."""
    def full(self, oracle_gain=0.2, cap=True):
        m = Decide().grid(oracle_gain)
        m['cap_reached'] = cap
        m['error'] = pd.Series(np.nan, index=m.index, dtype=object)   # empty string reads back from CSV as NaN
        m['total_transport'] = 1.0
        m['ligand_allocation_tv_mean'] = 0.2
        m['frac_ligand'] = 0.5
        return m

    def test_rules(self):
        m = self.full()
        with self.assertRaises(Exception):
            si.decide(m)                          # plan.md rule: capped baselines stop the summary
        dec, _ = si.decide(m, rule='checks')
        self.assertTrue(dec['R1_oracle_usable']); self.assertEqual(dec['runs_not_converged'], 0)
        self.assertEqual(dec['runs_at_cap'], len(m)); self.assertEqual(dec['convergence_rule'], 'checks')
        bad = (m.seed == 2) & (m.variant == 'q100') & (m.norm == 'raw') & (m.beta == 1.0)
        m.loc[bad, 'error'] = 'Budget exceeded: solver has not converged'
        m.loc[bad, ['pooled_spearman', 'total_transport']] = np.nan
        dec, _ = si.decide(m, rule='checks')
        self.assertFalse(dec['R1_cells']['raw|1']); self.assertTrue(dec['R1_cells']['raw|0.1'])
        self.assertEqual(dec['runs_not_converged'], 1)
        with self.assertRaises(Exception):
            si.decide(m, rule='nonsense')

    def test_stability(self):
        used = self.full()
        sub = lambda m: m[(m.seed == 1) & m.variant.isin(['none', 'q100', 'q100_perm'])].copy()
        r = si.stability_check(used, sub(self.full()))
        self.assertTrue(r['passed']); self.assertEqual(r['compared_solves'], 1 + 2 * 2 * 2)
        s = sub(self.full()); s['pooled_spearman'] += 5e-4
        self.assertTrue(si.stability_check(used, s)['passed'])
        s = sub(self.full()); s['total_transport'] *= 1 + 2e-3
        r = si.stability_check(used, s); self.assertFalse(r['passed'])
        self.assertAlmostEqual(r['max_abs_diff']['total_transport_relative'], 2e-3 / (1 + 2e-3))
        s = sub(self.full()); s.loc[s.index[2], 'error'] = 'Budget exceeded'; s.loc[s.index[2], 'total_transport'] = np.nan
        r = si.stability_check(used, s); self.assertFalse(r['passed']); self.assertEqual(r['failed_checks_long'], 1)
        s = sub(self.full()); s['seed'] = 9                     # no counterpart among the used runs
        with self.assertRaises(Exception):
            si.stability_check(used, s)


class Recheck(unittest.TestCase):
    def test_commands(self):
        import argparse, tempfile
        with tempfile.TemporaryDirectory() as d:
            d = Path(d); runs = d/'runs'
            lines = []
            for st, setting in (('main', '--eps 0.1 --rho 10'), ('sens', '--eps 0.2 --rho 1')):
                for name in ('s1_none', 's1_q100_raw', 's2_q100_raw'):
                    (runs/st/name).mkdir(parents=True); (runs/st/name/'metrics.csv').write_text('x')
                    lines.append(f'conda run -n inject-gpu python stage2_inject.py run --truth {runs}/truth_s1 {setting} '
                                 f'--nitermax 200000 --device cuda --output {runs}/{st}/{name} > {runs}/logs/{st}_{name}.log 2>&1')
            (runs/'commands.txt').write_text('\n'.join(lines) + '\n')
            ns = lambda **k: argparse.Namespace(**{**dict(commands=runs/'commands.txt', runs=str(runs), check_runs=str(d/'chk'),
                                                         nitermax=400000, output=d/'cmds.txt',
                                                         select=['s1_none', 's1_q100_raw'], require_all=True), **k})
            info = si.recheck_commands(ns())
            self.assertEqual((info['selected_lines'], info['original_nitermax'], info['device']), (2, 200000, 'cuda'))
            out = (d/'cmds.txt').read_text().splitlines()
            self.assertTrue(all('--nitermax 400000' in l and f'{d}/chk/main/' in l and f'> {d}/chk/logs/' in l
                                and f'{runs}/truth_s1' in l and '--eps 0.1 --rho 10' in l for l in out))
            with self.assertRaises(Exception):                       # cap must go up
                si.recheck_commands(ns(nitermax=100000, output=d/'c2.txt'))
            with self.assertRaises(Exception):                       # --runs spelled differently
                si.recheck_commands(ns(runs=str(runs) + '_x', output=d/'c3.txt'))
            with self.assertRaises(Exception):                       # a required stability run is absent
                si.recheck_commands(ns(select=['s1_none', 's1_uniform'], output=d/'c5.txt'))
            self.assertEqual(si.recheck_commands(ns(select=['s1_none', 's1_uniform'], require_all=False,
                                                    output=d/'c6.txt'))['selected'], ['s1_none'])
            (runs/'sens'/'s1_none'/'metrics.csv').unlink()
            with self.assertRaises(Exception):                       # an original run is missing
                si.recheck_commands(ns(output=d/'c4.txt'))


class Device(unittest.TestCase):
    def test_single_device(self):
        import pandas as pd
        import stage2_inject as si
        self.assertEqual(si.single_device(pd.DataFrame(dict(beta=[0.1]))), 'cpu')
        self.assertEqual(si.single_device(pd.DataFrame(dict(device=['cuda', 'cuda']))), 'cuda')
        with self.assertRaises(Exception):
            si.single_device(pd.DataFrame(dict(device=['cpu', 'cuda'])))


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
