"""Toy tests for isoform_lr_bench.py: ranks, metrics, jobs -> merge round trip, training overlap, verdict."""
import argparse
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

import isoform_lr_bench as ib


def toy_dataset():
    """Two tiers, two contrasts each; reference listed first in every contrast."""
    rows = []
    seqs = iter(['ACDEFGHIKL' * k for k in range(3, 40)])
    for tier in ib.TIERS:
        for c in range(2):
            partner = next(seqs)
            for i, role in enumerate(['reference', 'competitor', 'competitor', 'excluded_other_evidence']):
                v = next(seqs)
                rows.append(dict(pair_id=f'{tier}{c}{i}', tier=tier, contrast_id=f'{tier}:{c}', role_in_contrast=role,
                                 variant_id=f'G{c}-{i}', variant_gene=f'G{c}', variant_status='Displayed' if i == 1 else 'Described',
                                 partner_id=f'P{c}', pct_id_to_reference=100.0 if i == 0 else 60.0 + 10 * i,
                                 len_variant=len(v), len_partner=len(partner), n_tokens=len(v) + len(partner) + 3,
                                 seq_variant=v, seq_partner=partner))
    return pd.DataFrame(rows)


class Ranks(unittest.TestCase):
    def test_ties_share_average_rank(self):
        np.testing.assert_allclose(ib.ranks_desc([0.9, 0.5, 0.9, 0.1]), [1.5, 3, 1.5, 4])

    def test_contrast_metrics(self):
        d = toy_dataset()
        d = d[d.role_in_contrast.isin(ib.RANKED_ROLES)].copy()
        # contrast 0: reference top; contrast 1: reference last
        d['score'] = [0.9, 0.1, 0.2, 0.1, 0.5, 0.6] * 2
        t = ib.contrast_table(d[d.tier == 'H_isoform'])
        self.assertEqual(t.hit_at_1.tolist(), [1, 0])
        self.assertEqual(t.norm_rank.tolist(), [0.0, 1.0])
        self.assertAlmostEqual(t.delta.iloc[0], 0.7)

    def test_permutation_null_is_calibrated(self):
        t = pd.DataFrame(dict(k=[4] * 200, norm_rank=[0.5] * 200, hit_at_1=[0.25] * 200))
        p_nr, p_h1 = ib.permutation_p(t, 2000, np.random.default_rng(0))
        self.assertGreater(p_nr, 0.2); self.assertGreater(p_h1, 0.2)
        t['norm_rank'] = 0.0; t['hit_at_1'] = 1.0
        p_nr, p_h1 = ib.permutation_p(t, 2000, np.random.default_rng(0))
        self.assertLess(p_nr, 0.001); self.assertLess(p_h1, 0.001)


class RoundTrip(unittest.TestCase):
    def test_jobs_merge_evaluate(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            d = toy_dataset(); d.to_csv(tmp/'ds.csv', index=False)
            ib.jobs(argparse.Namespace(dataset=tmp/'ds.csv', output_dir=tmp/'jobs', limit_contrasts=0,
                                       include_longest=False, seed=0, exploratory=True))
            jobs = pd.read_csv(tmp/'jobs'/'plm_input.csv')
            self.assertEqual(len(jobs), 2 * len(d))
            # fake plm_score.py output: reference gets 0.8 in one orientation and 0.6 in the other
            idx = pd.read_csv(tmp/'jobs'/'index.csv')
            ref = set(d.loc[d.role_in_contrast == 'reference', 'pair_id'])
            sc = [0.8 if (p in ref and o == 'vp') else 0.6 if p in ref else 0.2 for p, o in zip(idx.pair_id, idx.orientation)]
            pd.DataFrame(dict(pair_id=idx.job_id, score=sc, n_tokens=0)).to_csv(tmp/'s.csv', index=False)
            (tmp/'s.csv.meta.json').write_text(json.dumps(dict(
                status='complete', input_sha256=ib.digest(tmp/'jobs'/'plm_input.csv'),
                checkpoint_sha256='x', scored_beyond_max_tokens=[])))
            ib.merge(argparse.Namespace(dataset=tmp/'ds.csv', jobs_dir=tmp/'jobs', scores=tmp/'s.csv',
                                        output=tmp/'m.csv', exploratory=True))
            m = pd.read_csv(tmp/'m.csv')
            self.assertTrue(np.allclose(m.loc[m.role_in_contrast == 'reference', 'score'], 0.7))
            with self.assertRaises(FileExistsError):
                ib.merge(argparse.Namespace(dataset=tmp/'ds.csv', jobs_dir=tmp/'jobs', scores=tmp/'s.csv',
                                            output=tmp/'m.csv', exploratory=True))
            summary, pc, pw, _ = ib.evaluate_frame(m, 500, 50, 0)
            self.assertEqual(summary['H_isoform']['hit_at_1'], 1.0)
            self.assertEqual(len(pc), 4)                      # excluded rows never ranked
            self.assertEqual(set(pw.win), {1.0})

    def test_merge_refuses_missing_jobs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            d = toy_dataset(); d.to_csv(tmp/'ds.csv', index=False)
            ib.jobs(argparse.Namespace(dataset=tmp/'ds.csv', output_dir=tmp/'jobs', limit_contrasts=0,
                                       include_longest=False, seed=0, exploratory=True))
            idx = pd.read_csv(tmp/'jobs'/'index.csv').iloc[:-1]
            pd.DataFrame(dict(pair_id=idx.job_id, score=0.5)).to_csv(tmp/'s.csv', index=False)
            (tmp/'s.csv.meta.json').write_text(json.dumps(dict(
                status='complete', input_sha256=ib.digest(tmp/'jobs'/'plm_input.csv'), scored_beyond_max_tokens=[])))
            with self.assertRaises(ValueError):
                ib.merge(argparse.Namespace(dataset=tmp/'ds.csv', jobs_dir=tmp/'jobs', scores=tmp/'s.csv',
                                            output=tmp/'m.csv', exploratory=True))

    def test_registered_dataset_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            toy_dataset().to_csv(Path(tmp)/'ds.csv', index=False)
            with self.assertRaises(ValueError):
                ib.load_dataset(Path(tmp)/'ds.csv', check_sha=True)


class Overlap(unittest.TestCase):
    def test_training_overlap_uses_exact_sequences(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            d = toy_dataset()
            first = d.iloc[0]
            f = tmp/ib_training_file()
            f.parent.mkdir(parents=True)
            pd.DataFrame(dict(query=[first.seq_partner], text=[first.seq_variant], label=['1'])).to_csv(f, index=False)
            rows = ib.training_overlap(d, tmp, 'humanV12')
            self.assertEqual(int(rows.train_pos.sum()), 1)
            self.assertTrue(bool(rows.train_pos.iloc[0]))


def ib_training_file():
    import lr_eval
    return lr_eval.TRAINING['humanV12'][0]


class Verdict(unittest.TestCase):
    def test_decision_rule(self):
        H = dict(p_hit_at_1=0.001, p_norm_rank=0.001, hit_at_1=0.6, mean_norm_rank=0.3)
        E = dict(hit_at_1=0.7, p_hit_at_1=0.001)
        base = dict(canonical=0.34, longest=0.38)
        self.assertEqual(ib.decide(H, E, dict(n_contrasts=40, p_norm_rank=0.01), base)['VERDICT_isoform_resolution'], 'PASS')
        self.assertEqual(ib.decide(H, E, dict(n_contrasts=40, p_norm_rank=0.2), base)['VERDICT_isoform_resolution'], 'FAIL')
        self.assertEqual(ib.decide(H, E, dict(n_contrasts=10, p_norm_rank=0.2), base)['VERDICT_isoform_resolution'], 'PASS')
        self.assertTrue(ib.decide(H, E, None, base)['VERDICT_isoform_resolution'].startswith('incomplete'))
        self.assertEqual(ib.decide(dict(H, hit_at_1=0.46), E, None, base)['C3_beats_trivial_baselines'], False)


if __name__ == '__main__':
    unittest.main()
