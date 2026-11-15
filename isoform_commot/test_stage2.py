import argparse
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

import plans
import stage2
from research import backend_sparse, counted, solve, solve_sparse, sparse_kkt_check, sparse_mass_check


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
            self.assertFalse(info['cap_reached'])
            self.assertLess(0, info['niter'])
            self.assertLessEqual(info['niter'], info['nitermax'])
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

    def test_kkt_check_large_index_space(self):
        # whole-slide sizes: row * n_cols exceeds int32 (csgraph returns int32 node indices)
        rng = np.random.default_rng(9)
        C = rng.uniform(0, 1, (6, 5))
        a, b = rng.uniform(.1, 1, 6), rng.uniform(.1, 1, 5)
        r, c = np.nonzero(np.isfinite(C))
        T, info = solve_sparse(r, c, C[r, c], a, b, 0.1, 10.0)
        s, R, K, off = max(a.sum(), b.sum()), 60000, 60000, 59990
        aa, bb = np.zeros(R), np.zeros(K)
        aa[off:off+6], bb[off:off+5] = a/s, b/s
        out = sparse_kkt_check(r + off, c + off, C[r, c], T/s, aa, bb, 0.1, 10.0)
        self.assertLessEqual(out['kkt_tol'], 1e-9)

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
        early, iterations = counted(backend_sparse(), a/s, b/s, coo, 0.1, 10.0, nitermax=3, stopthr=1e-11)
        self.assertEqual(iterations, [4])                           # cap reached: nitermax + 1
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
        module = commot_baseline.official_cot()
        combine, collective = module.cot_combine_sparse, module.cot_sparse
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
        P3, iterations = commot_baseline.solve_counted(module, 'cot_sparse', S, D, A, M, cut, eps_p=0.1, rho=10.0,
                                                       nitermax=50000, stopthr=1e-11)
        self.assertEqual(len(iterations), 1)                        # one collective solve, counted
        self.assertIs(module.unot.__name__, 'unot')                 # original restored
        for key in P1:
            np.testing.assert_array_equal(P1[key].toarray(), P3[key].toarray())   # counting changes nothing
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
            head = subprocess.run(['git', '-C', str(Path(stage2.__file__).parent), 'rev-parse', 'HEAD'],
                                  capture_output=True, text=True, check=True).stdout.strip()
            self.assertEqual(manifest['launch']['git_head'], head)
            self.assertIsInstance(manifest['launch']['git_uncommitted_changes'], bool)
            self.assertLessEqual(plan['D'].max(), 250.0)
            self.assertGreater(summary.total_T_contact.iloc[0], 0)
            self.assertEqual(summary.check_rows_checked.iloc[0], summary.check_rows.iloc[0])


class Report(unittest.TestCase):
    def make(self, root, partial, input_dir, niter=100):
        cut = stage2.COMMOT_TUTORIAL_CUTOFF_UM
        ours, commot, ev = root/'A', root/'B', root/'eval'
        for folder in (ours, commot, ev):
            folder.mkdir()
        (ours/'manifest.json').write_text(json.dumps(dict(method='research sparse (stage 2)',
            config=dict(nitermax=50000, stopthr=1e-11))))
        pd.DataFrame(dict(case=['case_000', 'case_001', 'case_002'], r0_um=[50.0, 100.0, 200.0], eps=0.1,
                          rho=10.0, diffusion_cutoff_um=cut, frac_ligand=.5, frac_receptor=.4, total_T=1.0,
                          mean_distance=10.0, solve_seconds=1.0, niter=[niter, 50001, niter], nitermax=50000,
                          cap_reached=[False, True, False], check_kkt_tol=1e-11)).to_csv(ours/'summary.csv', index=False)
        (commot/'manifest.json').write_text(json.dumps(dict(method='COMMOT original')))
        pd.DataFrame(dict(case=['case_000'], r0=[cut], dis_thr_um=[cut], eps=0.1, rho=10.0, scheme='collective',
                          nitermax=50000, stopthr=1e-11, frac_ligand=.5, frac_receptor=.4, total_T=1.0,
                          mean_distance=10.0, solve_seconds=1.0, check_passed=True, check_kkt_tol=1e-11, niter=niter,
                          solves=1, solves_at_cap=0, cap_reached=False)).to_csv(commot/'summary.csv', index=False)
        pd.DataFrame(dict(receptor=['R1', 'R2', 'R3'], top_tf=10, status='evaluated')).to_csv(ev/'target_sets.csv',
                                                                                              index=False)
        common = dict(shape=(4, 4), a=np.ones(4), b=np.ones(4), spots=['s0', 's1'], ligands=['L1', 'L2'],
                      receptors=['R1', 'R2'], pair_type=np.array(['diffusion', 'contact']))
        for k, far in enumerate([False, True, False]):            # case_001 carries mass beyond the cutoff
            d = [10.0, cut + 50.0 if far else 20.0, 100.0]
            plans.save_sparse(ours/f'case_{k:03d}.npz', [0, 1, 2], [0, 1, 2], [0.1, 0.2, 0.0], [1.0, 0.5, 0.2],
                              common['shape'], common['a'], common['b'], common['spots'], common['ligands'],
                              common['receptors'], D=np.array(d), pair=np.array([0, 0, 1]), pair_type=common['pair_type'])
        plans.save_sparse(commot/'case_000.npz', [0, 1], [0, 1], [10.0, cut], [1.0, 0.5], common['shape'], common['a'],
                          common['b'], common['spots'], common['ligands'], common['receptors'])
        rows = [dict(results=str(root/f), case=c, level='receptor', top_tf=10, target_score='mean',
                     score_mode='no_ppi', receptor=r, partial_spearman=v, spearman=v, receptor_only_spearman=0.0)
                for (f, c), values in partial.items() for r, v in zip(['R1', 'R2', 'R3'], values)]
        pd.DataFrame(rows).to_csv(ev/'case_metrics.csv', index=False)
        (ev/'eval_manifest.json').write_text(json.dumps(dict(status='complete', results=[str(ours), str(commot)],
            primary_top_tf=10, target_score=['mean'], input=str(input_dir))))
        return ev

    def test_r0_rule_paired_sign_test_classes_and_confirmation(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root/'in').mkdir()
            pd.DataFrame(dict(ligand=['L1', 'L2', 'L1', 'L2'], receptor=['R1', 'R2', 'R3', 'R3'],
                              type=['diffusion', 'contact', 'diffusion', 'contact'])).to_csv(root/'in'/'lr.csv', index=False)
            # R1, R2 defined at every r0: means 0.3, 0.3, 0.2 -> tie between 50 and 100, rule picks 50.
            # R3 undefined at r0 = 200: excluded from the selection although it favours 100.
            post = {('A', 'case_000'): [.2, .4, .0], ('A', 'case_001'): [.4, .2, .9],
                    ('A', 'case_002'): [.2, .2, np.nan], ('B', 'case_000'): [.1, .1, .1]}
            (root/'post').mkdir()
            ev = self.make(root/'post', post, root/'in')
            stage2.report(argparse.Namespace(eval=str(ev), output=str(root/'post'/'report'), posterior_report=None))
            out = json.loads((root/'post'/'report'/'report.json').read_text())
            self.assertIn('git_head', out['launch'])
            self.assertEqual((out['role'], out['selected_r0_um'], out['receptors_for_r0_selection']), ('selection', 50.0, 2))
            self.assertFalse(out['all_converged'])
            conv = pd.read_csv(root/'post'/'report'/'convergence.csv')
            self.assertEqual(list(conv.converged), [True, False, True, True])   # mass beyond the cutoff
            self.assertEqual(list(conv.forbidden_mass), [0.0, 0.5, 0.0, 0.0])
            self.assertTrue(conv.cap_reached.iloc[1])                    # recorded; not by itself a failure
            paired = pd.read_csv(root/'post'/'report'/'paired_receptors.csv').set_index('receptor_class')
            self.assertEqual(tuple(paired.loc['all', ['receptors', 'ours_higher', 'baseline_higher']]), (3, 2, 1))
            self.assertAlmostEqual(paired.loc['all', 'sign_test_p'], 1.0)
            self.assertEqual(paired.loc['Secreted', 'receptors'], 1)    # R1 diffusion only
            self.assertEqual(paired.loc['Contact', 'receptors'], 1)     # R2 contact only
            self.assertEqual(paired.loc['both', 'receptors'], 1)        # R3 both types
            pos = pd.read_csv(root/'post'/'report'/'positive_fraction.csv')
            sel = pos[pos.selected_r0].iloc[0]
            self.assertEqual((sel.r0_um, sel.partial_positive, sel.receptors_defined), (50.0, 2, 3))
            # anterior: the same rule picks 100 here, so the posterior r0 (50) is not confirmed
            ant = {('A', 'case_000'): [.1, .1, .1], ('A', 'case_001'): [.3, .3, .3],
                   ('A', 'case_002'): [.2, .2, .2], ('B', 'case_000'): [.0, .0, .0]}
            (root/'ant').mkdir()
            ev = self.make(root/'ant', ant, root/'in')
            stage2.report(argparse.Namespace(eval=str(ev), output=str(root/'ant'/'report'),
                                             posterior_report=str(root/'post'/'report'/'report.json')))
            out = json.loads((root/'ant'/'report'/'report.json').read_text())
            c = out['confirmation']
            self.assertEqual((out['role'], out['selected_r0_um']), ('confirmation', 50.0))
            self.assertEqual((c['rule_selected_r0_here'], c['confirmed'], c['rank_of_posterior_r0']), (100.0, False, 3))
            self.assertAlmostEqual(c['gap_to_best'], 0.2)
            self.assertFalse(c['strict_confirmed'])
            self.assertAlmostEqual(c['delta'], 0.2)
            self.assertAlmostEqual(c['bootstrap_se'], 0.0)             # identical receptors: no spread
            self.assertFalse(c['tolerance_confirmed'])

    def test_tolerance_confirmation_bootstrap(self):
        a = pd.DataFrame({50.0: [.1, .5, .3, .2], 100.0: [.2, .4, .3, .3]})
        out = stage2.tolerance_confirmation(a, 100.0, 50.0)
        self.assertAlmostEqual(out['delta'], 0.025)
        self.assertGreater(out['bootstrap_se'], out['delta'])          # within one standard error
        self.assertTrue(out['tolerance_confirmed'])
        self.assertEqual(out, stage2.tolerance_confirmation(a, 100.0, 50.0))   # fixed seed


class Stability(unittest.TestCase):
    def test_50k_vs_200k_and_determinism(self):
        sp = hex_slide(5, 10)
        rng = np.random.default_rng(12)
        x = pd.DataFrame(rng.poisson(3, (len(sp), 3)).astype(float) + 1, index=sp.index, columns=['L1', 'R1', 'R2'])
        lr = pd.DataFrame(dict(ligand=['L1', 'L1'], receptor=['R1', 'R2'], type=['diffusion', 'contact']))
        cfg = dict(stage2.DEFAULTS, r0_um=[50.0], coordinate_unit='micrometer', expression_level='gene',
                   expression_state='processed')
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            write_input(root/'in', sp, x, lr, cfg)
            for name, n in [('k50', 50000), ('k200', 200000), ('again', 50000)]:
                stage2.run(argparse.Namespace(input=str(root/'in'), output=str(root/name), r0_um=None, eps=0.1, rho=10.0,
                                              diffusion_cutoff_um=None, nitermax=n, stopthr=None))
            args = dict(at_50k=[str(root/'k50')], at_200k=[str(root/'k200')], same_a=[str(root/'k50')],
                        same_b=[str(root/'again')])
            stage2.stability(argparse.Namespace(output=str(root/'s1'), **args))
            stage2.convergence(argparse.Namespace(results=[str(root/'k50')], output=str(root/'conv')))
            conv = pd.read_csv(root/'conv'/'convergence.csv').iloc[0]
            self.assertTrue(conv.converged and conv.forbidden_mass == 0 and not conv.cap_reached)
            v = pd.read_csv(root/'s1'/'stability.csv').iloc[0]
            self.assertTrue(v.stable and v.rel_total_difference == 0 and v.receptors == 2)
            self.assertTrue(pd.read_csv(root/'s1'/'determinism.csv').identical.all())
            with np.load(root/'k200'/'case_000.npz', allow_pickle=False) as z:   # 1% more transport: fails
                data = {k: z[k] for k in z.files}
            data['T'] = data['T'] * 1.01
            np.savez_compressed(root/'k200'/'case_000.npz', **data)
            stage2.stability(argparse.Namespace(output=str(root/'s2'), **args))
            v = pd.read_csv(root/'s2'/'stability.csv').iloc[0]
            self.assertFalse(v.stable)
            self.assertAlmostEqual(v.rel_total_difference, 0.01 / 1.01)
            self.assertFalse(json.loads((root/'s2'/'stability.json').read_text())['all_stable'])


class CommotRule(unittest.TestCase):
    def test_commot_target_sets(self):
        from evaluate_targets import commot_target_sets
        rt = pd.DataFrame(dict(receptor=['Fgfr1', 'Fgfr1', 'Acvr2a,Acvr1', 'Tgfbr2'], pathway='p',
                               tf=['T1', 'T2', 'T3', 'T4'], weight=[1e-3, 1e-6, 1e-4, 1e-3],
                               source=['REACTOME', 'REACTOME', 'REACTOME', 'KEGG']))
        tt = pd.DataFrame(dict(tf=['T1', 'T1', 'T2', 'T3', 'T4'], target=['G1', 'Fgfr1', 'G2', 'G3', 'G4']))
        hvg = {'G1': True, 'Fgfr1': True, 'G2': True, 'G3': True, 'G4': True}
        symbol = {s: s for s in ['Fgfr1', 'Acvr1', 'Acvr2a', 'Tgfbr2']}
        sets = commot_target_sets(['Fgfr1', 'ACVR', 'Tgfbr2'], symbol, rt, tt, hvg,
                                  subunits={'ACVR': ['Acvr1', 'Acvr2a']}).set_index('receptor')
        self.assertEqual(sets.loc['Fgfr1', 'tfs'], 'T1')                 # T2 below PPR 1e-5
        self.assertEqual(sets.loc['Fgfr1', 'targets'], 'Fgfr1;G1')       # LR genes kept
        self.assertEqual(sets.loc['ACVR', 'targets'], 'G3')              # complex matched as a subunit set
        self.assertFalse(sets.loc['ACVR', 'ordered_match'])              # order differs from scSeqComm
        self.assertEqual(sets.loc['Tgfbr2', 'status'], 'no_receptor_tf_link')   # KEGG only
        sets = commot_target_sets(['Fgfr1'], symbol, rt, tt, {'G1': False, 'Fgfr1': False})
        self.assertEqual(sets.status.iloc[0], 'no_hvg_target')


class Robustness(unittest.TestCase):
    def test_subsample_and_robustness(self):
        sp = hex_slide(6, 12)
        rng = np.random.default_rng(11)
        x = pd.DataFrame(rng.poisson(3, (len(sp), 3)).astype(float) + 1, index=sp.index, columns=['L1', 'R1', 'R2'])
        lr = pd.DataFrame(dict(ligand=['L1', 'L1'], receptor=['R1', 'R2'], type=['diffusion', 'contact'],
                               pathway_name='p'))
        cfg = dict(stage2.DEFAULTS, r0_um=[100.0], diffusion_cutoff_um=250.0, coordinate_unit='micrometer',
                   expression_level='gene', expression_state='processed')
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            write_input(root/'in', sp, x, lr, cfg)
            lr.assign(reason='selected').to_csv(root/'in'/'filter_audit.csv', index=False)
            prep = json.loads((root/'in'/'preparation.json').read_text())
            prep['prepared_sha256']['filter_audit.csv'] = stage2.digest(root/'in'/'filter_audit.csv')
            (root/'in'/'preparation.json').write_text(json.dumps(prep))
            stage2.subsample(argparse.Namespace(input=str(root/'in'), output=str(root/'sub'), fraction=0.5, seed=3))
            again = root/'sub2'
            stage2.subsample(argparse.Namespace(input=str(root/'in'), output=str(again), fraction=0.5, seed=3))
            self.assertEqual((root/'sub'/'spatial.csv').read_text(), (again/'spatial.csv').read_text())  # seeded
            kept = pd.read_csv(root/'sub'/'spatial.csv', index_col=0)
            self.assertEqual(len(kept), round(0.5*len(sp)))
            for name, src in [('full', 'in'), ('subrun', 'sub')]:
                stage2.run(argparse.Namespace(input=str(root/src), output=str(root/name), r0_um=None, eps=0.1, rho=10.0,
                                              diffusion_cutoff_um=stage2.COMMOT_TUTORIAL_CUTOFF_UM))
            stage2.robustness(argparse.Namespace(full=[str(root/'full')], subsampled=[str(root/'subrun')],
                                                 output=str(root/'rob')))
            table = pd.read_csv(root/'rob'/'robustness_receptors.csv')
            self.assertEqual(set(table.receptor), {'R1', 'R2'})
            self.assertTrue(((table.fraction == 0.5) & (table.seed == 3)).all())
            self.assertTrue(table.spearman.between(-1, 1).all())


class Validation(unittest.TestCase):
    """Plan section 9: flat and inverse kernels, coordinate shuffle."""

    def test_diffusion_cost_kernels(self):
        d = np.array([0.0, 100.0, 200.0, 300.0])
        np.testing.assert_allclose(stage2.diffusion_cost(d, 'r0', r0=100.0), np.log1p(d / 100.0))
        np.testing.assert_array_equal(stage2.diffusion_cost(d, 'flat'), np.zeros(4))
        np.testing.assert_allclose(stage2.diffusion_cost(d, 'inverse', d_ref=100.0), [0, 0, np.log(2), np.log(3)])
        with self.assertRaises(ValueError):
            stage2.diffusion_cost(d, 'inverse')
        sp = hex_slide(6, 12)
        self.assertAlmostEqual(stage2.min_spot_distance(sp), 100.0, places=6)

    def test_flat_inverse_runs_and_shuffle(self):
        sp = hex_slide(6, 12)
        rng = np.random.default_rng(5)
        x = pd.DataFrame(rng.poisson(3, (len(sp), 3)).astype(float) + 1, index=sp.index, columns=['L1', 'R1', 'R2'])
        lr = pd.DataFrame(dict(ligand=['L1', 'L1'], receptor=['R1', 'R2'], type=['diffusion', 'contact'],
                               pathway_name='p'))
        cfg = dict(stage2.DEFAULTS, r0_um=[100.0, 200.0], diffusion_cutoff_um=250.0, coordinate_unit='micrometer',
                   expression_level='gene', expression_state='processed')
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            write_input(root/'in', sp, x, lr, cfg)
            lr.assign(reason='selected').to_csv(root/'in'/'filter_audit.csv', index=False)
            prep = json.loads((root/'in'/'preparation.json').read_text())
            prep['prepared_sha256']['filter_audit.csv'] = stage2.digest(root/'in'/'filter_audit.csv')
            (root/'in'/'preparation.json').write_text(json.dumps(prep))
            base = dict(input=str(root/'in'), eps=0.2, rho=1.0, diffusion_cutoff_um=250.0)
            stage2.run(argparse.Namespace(output=str(root/'A'), r0_um=[200.0], kernel='r0', **base))
            stage2.run(argparse.Namespace(output=str(root/'F'), r0_um=None, kernel='flat', **base))
            stage2.run(argparse.Namespace(output=str(root/'I'), r0_um=None, kernel='inverse', **base))
            with self.assertRaises(ValueError):
                stage2.run(argparse.Namespace(output=str(root/'bad'), r0_um=[100.0], kernel='flat', **base))
            pa, pf, pi = (plans.load(root/k/'case_000.npz') for k in 'AFI')
            for k in ['row', 'col', 'D']:                                  # same connections, only the cost differs
                np.testing.assert_array_equal(pa[k], pf[k]); np.testing.assert_array_equal(pa[k], pi[k])
            diffusion = pa['pair_type'][pa['pair']] == 'diffusion'
            self.assertTrue((pf['C'] == 0).all())
            np.testing.assert_allclose(pa['C'][diffusion], np.log1p(pa['D'][diffusion] / 200.0))
            np.testing.assert_allclose(pi['C'][diffusion], np.maximum(np.log(np.maximum(pi['D'][diffusion], 100.0) / 100.0), 0),
                                       atol=1e-9)
            self.assertTrue((pi['C'][~diffusion] == 0).all() and (pa['C'][~diffusion] == 0).all())
            for k in 'FI':
                summary = pd.read_csv(root/k/'summary.csv')
                self.assertEqual(len(summary), 1)
                self.assertTrue(summary.r0.isna().all())
                manifest = json.loads((root/k/'manifest.json').read_text())
                self.assertEqual(manifest['kernel'], {'F': 'flat', 'I': 'inverse'}[k])
            self.assertAlmostEqual(json.loads((root/'I'/'manifest.json').read_text())['d_ref_um'], 100.0, places=6)
            self.assertEqual(pd.read_csv(root/'F'/'summary.csv').family.tolist(), ['flat'])

            stage2.shuffle(argparse.Namespace(input=str(root/'in'), output=str(root/'s1'), seed=1))
            stage2.shuffle(argparse.Namespace(input=str(root/'in'), output=str(root/'s1b'), seed=1))
            stage2.shuffle(argparse.Namespace(input=str(root/'in'), output=str(root/'s2'), seed=2))
            self.assertEqual((root/'s1'/'spatial.csv').read_text(), (root/'s1b'/'spatial.csv').read_text())
            self.assertNotEqual((root/'s1'/'spatial.csv').read_text(), (root/'s2'/'spatial.csv').read_text())
            self.assertEqual((root/'s1'/'expression.csv').read_bytes(), (root/'in'/'expression.csv').read_bytes())
            moved = pd.read_csv(root/'s1'/'spatial.csv', index_col=0)
            original = pd.read_csv(root/'in'/'spatial.csv', index_col=0)
            self.assertEqual(list(moved.index), list(original.index))       # spots keep their order and expression
            perm = np.random.default_rng(1).permutation(len(original))
            np.testing.assert_array_equal(moved.to_numpy(), original.to_numpy()[perm])   # rows move as a whole
            self.assertEqual(list(moved.dtypes), list(original.dtypes))
            stage2.run(argparse.Namespace(output=str(root/'P'), r0_um=[200.0], kernel='r0',
                                          **dict(base, input=str(root/'s1'))))
            pp = plans.load(root/'P'/'case_000.npz')
            self.assertEqual(list(pp['spots']), list(pa['spots']))
            np.testing.assert_array_equal(pp['a'], pa['a'])                    # budgets unchanged
            self.assertFalse(np.array_equal(np.sort(pp['row']*10**6 + pp['col']), np.sort(pa['row']*10**6 + pa['col'])))
            self.assertEqual(json.loads((root/'P'/'manifest.json').read_text())['shuffle']['seed'], 1)
            with self.assertRaises(ValueError):
                stage2.shuffle(argparse.Namespace(input=str(root/'s1'), output=str(root/'s1x'), seed=1))


class R0Validation(unittest.TestCase):
    """Plan section 9: condition codes, evaluation of shuffled results, criteria (a)-(d), negative control."""

    def table(self, a, f, i, b, n=40, seed=0):
        rng = np.random.default_rng(seed)
        noise = lambda: rng.normal(0, 0.01, n)
        idx = [f'R{k}' for k in range(n)]
        return pd.DataFrame({'A_50': 0.1 + noise(), 'A_200': a + noise(), 'F': f + noise(), 'I': i + noise(),
                             'B': b + noise(), 'P1': 0.0 + noise(), 'P2': 0.0 + noise(), 'P3': 0.0 + noise()}, index=idx)

    def test_bootstrap_is_seeded_and_criteria_verdicts(self):
        import r0_validation as rv
        x = np.random.default_rng(3).normal(0.1, 0.05, 30)
        self.assertEqual(rv.bootstrap_mean(x), rv.bootstrap_mean(x))
        b = rv.bootstrap_mean(x)
        self.assertTrue(b['ci_low'] < b['mean'] < b['ci_high'])
        classes = {f'R{k}': 'Secreted' if k < 20 else 'Contact' if k < 35 else 'both' for k in range(40)}
        t = self.table(0.30, 0.20, 0.20, 0.20)
        t.loc[t.index[20:35], 'A_200'] = t.loc[t.index[20:35], 'F'] + 0.001        # contact receptors: no gain
        summary, receptors, common = rv.difference_tables(t, 200.0, classes)
        c = rv.criteria(summary, receptors)
        self.assertEqual([c[k]['met'] for k in 'abcd'], [True, True, True, True])
        self.assertEqual(c['verdict'], 'useful')
        self.assertEqual((c['d']['secreted_n'], c['d']['contact_n'], c['d']['both_n']), (20, 15, 5))
        self.assertGreater(c['d']['secreted_mean'], c['d']['contact_mean'])
        t = self.table(0.20, 0.30, 0.30, 0.30)                                       # A worse everywhere
        c = rv.criteria(*rv.difference_tables(t, 200.0, classes)[:2])
        self.assertFalse(c['a']['met'] or c['b']['met'] or c['c']['met'])
        t = self.table(0.30, 0.20, 0.20, 0.40)                                       # better than F and I, worse than B
        c = rv.criteria(*rv.difference_tables(t, 200.0, classes)[:2])
        self.assertTrue(c['a']['met'] and c['c']['met'] and not c['b']['met'])
        self.assertTrue(c['b']['I_met'] and not c['b']['B_met'])
        self.assertTrue(c['verdict'].startswith('partly met: a'))
        t = self.table(0.20, 0.20, 0.20, 0.20)                                       # interval includes 0
        c = rv.criteria(*rv.difference_tables(t, 200.0, classes)[:2])
        self.assertFalse(c['a']['met'])

    def test_negative_control_and_curve(self):
        import r0_validation as rv
        t = self.table(0.30, 0.2, 0.2, 0.2)
        table, verdict = rv.negative_control(t)
        self.assertFalse(verdict['insensitive'])                                     # A(200) clearly above shuffled
        self.assertAlmostEqual(verdict['mean_difference'], 0.30, places=2)
        t[['P1', 'P2', 'P3']] = t[['A_200']].to_numpy() + np.random.default_rng(1).normal(0, .01, (40, 3))
        self.assertTrue(rv.negative_control(t)[1]['insensitive'])
        self.assertEqual(rv.negative_control(t.drop(columns='P3')), (None, None))
        t.loc['R0', 'F'] = np.nan                                                    # R0 leaves the common set
        curve, common = rv.curve_table(t, [50.0, 200.0], pd.Series(0.05, index=t.index))
        self.assertEqual(len(common), 39)
        self.assertEqual(curve.condition.tolist(), ['A', 'A', 'F', 'I', 'B', 'receptor_only'])
        self.assertTrue((curve.n == 39).all())

    def test_condition_codes_and_report_on_real_runs(self):
        import commot_baseline
        import evaluate_targets
        import r0_validation as rv
        sp = hex_slide(5, 10)
        rng = np.random.default_rng(7)
        x = pd.DataFrame(rng.poisson(3, (len(sp), 4)).astype(float) + 1, index=sp.index, columns=['L1', 'R1', 'R2', 'R3'])
        lr = pd.DataFrame(dict(ligand=['L1', 'L1', 'L1'], receptor=['R1', 'R2', 'R3'],
                               type=['diffusion', 'contact', 'diffusion'], pathway_name='p'))
        cfg = dict(stage2.DEFAULTS, coordinate_unit='micrometer', expression_level='gene', expression_state='processed')
        cut = stage2.COMMOT_TUTORIAL_CUTOFF_UM
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            write_input(root/'in', sp, x, lr, cfg)
            lr.assign(reason='selected').to_csv(root/'in'/'filter_audit.csv', index=False)
            prep = json.loads((root/'in'/'preparation.json').read_text())
            prep['prepared_sha256']['filter_audit.csv'] = stage2.digest(root/'in'/'filter_audit.csv')
            (root/'in'/'preparation.json').write_text(json.dumps(prep))
            main = dict(eps=0.1, rho=10.0, diffusion_cutoff_um=cut)
            run = lambda name, src='in', **kw: stage2.run(argparse.Namespace(input=str(root/src), output=str(root/name),
                                                                             **dict(main, **kw)))
            run('A', r0_um=[100.0, 200.0], kernel='r0')
            run('F', r0_um=None, kernel='flat')
            run('I', r0_um=None, kernel='inverse')
            for seed in (1, 2, 3):
                stage2.shuffle(argparse.Namespace(input=str(root/'in'), output=str(root/f's{seed}'), seed=seed))
                run(f'P{seed}', src=f's{seed}', r0_um=[200.0], kernel='r0')
            commot_baseline.run(root/'in', root/'B', dis_thr_um=[cut], eps=0.1, rho=10.0, weights=(1, 0, 0, 0),
                                nitermax=50000, stopthr=1e-11)
            folders = ['A', 'F', 'I', 'P1', 'P2', 'P3', 'B']
            codes = {f: stage2.condition_cases(root/f).condition.tolist() for f in folders}
            self.assertEqual(codes, dict(A=['A', 'A'], F=['F'], I=['I'], P1=['P1'], P2=['P2'], P3=['P3'], B=['B']))
            for f in folders:                                    # shuffled results evaluate against the unshuffled input
                evaluate_targets.load_results(root/f, root/'in')
            with self.assertRaises(ValueError):
                evaluate_targets.load_results(root/'A', root/'s1')
            for f in ('F', 'P1'):
                case = stage2.condition_cases(root/f).to_dict('records')[0]
                case['forbidden_mass'] = stage2.forbidden_mass(case)
                self.assertTrue(stage2.converged(case))
            ev = root/'eval'
            ev.mkdir()
            receptors = ['R1', 'R2', 'R3']
            pd.DataFrame(dict(receptor=receptors, top_tf=10, status='evaluated')).to_csv(ev/'target_sets.csv', index=False)
            value = {('A', 'case_000'): [.2, .1, .3], ('A', 'case_001'): [.4, .1, .5], ('F', 'case_000'): [.1, .1, .2],
                     ('I', 'case_000'): [.2, .1, .2], ('B', 'case_000'): [.3, .1, .1], ('P1', 'case_000'): [.0, .1, .0],
                     ('P2', 'case_000'): [.1, .0, .0], ('P3', 'case_000'): [.0, .0, .1]}
            rows = [dict(results=str(root/f), case=c, level='receptor', top_tf=10, target_score='mean', score_mode='no_ppi',
                         receptor=r, partial_spearman=v, spearman=v, receptor_only_spearman=0.05)
                    for (f, c), values in value.items() for r, v in zip(receptors, values)]
            pd.DataFrame(rows).to_csv(ev/'case_metrics.csv', index=False)
            (ev/'eval_manifest.json').write_text(json.dumps(dict(status='complete', results=[str(root/f) for f in folders],
                primary_top_tf=10, target_score=['mean'], input=str(root/'in'))))
            rv.run(argparse.Namespace(eval=str(ev), output=str(root/'val_post'), posterior=None))
            post = json.loads((root/'val_post'/'validation.json').read_text())
            self.assertEqual((post['role'], post['r0_star_um'], post['receptors_curve']), ('selection', 200.0, 3))
            self.assertEqual(post['receptor_classes'], dict(Secreted=2, Contact=1))
            self.assertAlmostEqual(post['criteria']['a']['mean_difference'], (0.3 + 0.0 + 0.3) / 3)
            self.assertAlmostEqual(post['criteria']['d']['secreted_mean'], 0.3)
            self.assertAlmostEqual(post['criteria']['d']['contact_mean'], 0.0)
            self.assertAlmostEqual(post['negative_control']['mean_difference'], (0.4 + 0.1 + 0.5) / 3 - 0.1 / 3 * 3 / 3, places=6)
            for name in ['curve.csv', 'differences.csv', 'differences_receptors.csv', 'negative_control.csv', 'curve.png',
                         'differences.png', 'r0_selection.csv']:
                self.assertTrue((root/'val_post'/name).exists(), name)
            rv.run(argparse.Namespace(eval=str(ev), output=str(root/'val_ant'), posterior=str(root/'val_post'/'validation.json')))
            ant = json.loads((root/'val_ant'/'validation.json').read_text())
            self.assertEqual((ant['role'], ant['r0_star_um'], ant['criteria_role']), ('confirmation', 200.0, 'confirmatory (plan 9.5)'))
            stage2.report(argparse.Namespace(eval=str(ev), output=str(root/'rep'), posterior_report=None))   # section 8
            self.assertEqual(set(pd.read_csv(root/'rep'/'conditions.csv').condition), {'A', 'B'})          # ignores F, I, P


class Pairs(unittest.TestCase):
    def test_each_slide_keeps_its_own_pairs(self):
        sp = hex_slide(4, 8)
        x = pd.DataFrame(1.0, index=sp.index, columns=['L1', 'L2', 'R1'])
        cfg = dict(stage2.DEFAULTS, coordinate_unit='micrometer', expression_level='gene', expression_state='processed')
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for name, keep in [('p', [0, 1]), ('a', [1, 2])]:
                lr = pd.DataFrame(dict(ligand=['L1', 'L2', 'L1'], receptor=['R1', 'R1', 'R2'],
                                       type=['diffusion', 'contact', 'contact'], pathway_name='x'))
                write_input(root/name, sp, x, lr.iloc[keep], cfg)
                audit = lr.assign(reason=['selected' if i in keep else 'subunit_below_min_cell_pct' for i in range(3)])
                audit.to_csv(root/name/'filter_audit.csv', index=False)
            stage2.pairs(argparse.Namespace(posterior=str(root/'p'), anterior=str(root/'a'), output=str(root/'out')))
            counts = pd.read_csv(root/'out'/'counts.csv').set_index('type')
            self.assertEqual(tuple(counts.loc['all', ['posterior', 'anterior', 'intersection', 'union']]), (2, 2, 1, 3))
            union = pd.read_csv(root/'out'/'pairs_union.csv')
            only_p = union[union.in_posterior & ~union.in_anterior].iloc[0]
            self.assertEqual((only_p.ligand, only_p.anterior_filter), ('L1', 'subunit_below_min_cell_pct'))
            edges = pd.read_csv(root/'out'/'spot_pairs.csv')
            self.assertTrue((edges.max_distance_within <= edges.cutoff_um).all())


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
