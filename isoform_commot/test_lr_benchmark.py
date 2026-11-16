"""Toy tests for lr_benchmark.py (negatives) and lr_eval.py (metrics, ranking, consistency)."""
from collections import Counter
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

import lr_benchmark as lb
import lr_eval as le


class ToyEntities:
    """Single genes only; E = positives plus one extra annotated pair."""

    def __init__(self, positives, extra=()):
        self.e = {frozenset(p) for p in positives} | {frozenset(p) for p in extra}

    def sources(self, a, b):
        return {'toy'} if frozenset((a, b)) in self.e else set()

    def excluded(self, a, b):
        return bool(self.sources(a, b))

    def has_distinct_pair(self, a, b):
        return a != b

    def accessions(self, x):
        return [f'acc_{x}']

    def lengths(self, x):
        return [100]


def toy_positives(n_lig=12, n_rec=10, per=2, seed=0):
    rng = np.random.default_rng(seed)
    pairs = set()
    for i in range(n_lig):
        for j in rng.choice(n_rec, size=per, replace=False):
            pairs.add((f'L{i}', f'R{j}'))
    pos = pd.DataFrame(sorted(pairs), columns=['ligand', 'receptor'])
    pos['classes'] = 'Secreted'
    return pos


class Negatives(unittest.TestCase):
    def test_degree_preserving(self):
        pos = toy_positives()
        ent = ToyEntities(zip(pos.ligand, pos.receptor), extra=[('L0', 'R9')])
        for seed in (1, 2, 3):
            kept, stats = lb.degree_negatives(pos, ent, seed)
            self.assertEqual(stats['dropped_still_excluded'], 0)
            self.assertEqual(Counter(l for l, _ in kept), Counter(pos.ligand))
            self.assertEqual(Counter(r for _, r in kept), Counter(pos.receptor))
            self.assertFalse(any(ent.excluded(l, r) for l, r in kept))
            self.assertEqual(len({frozenset(e) for e in kept}), len(kept))
            self.assertGreaterEqual(stats['successful_swaps'], lb.MIN_SWAPS * len(pos))
        a, _ = lb.degree_negatives(pos, ent, 1)
        b, _ = lb.degree_negatives(pos, ent, 1)
        self.assertEqual(a, b)

    def test_candidates(self):
        pos = toy_positives()
        ent = ToyEntities(zip(pos.ligand, pos.receptor), extra=[('L0', 'R9')])
        cands, ligands, receptors, removed, annotated = lb.candidates(pos, ent)
        self.assertEqual(len(cands) + len(annotated), len(ligands) * len(receptors))
        self.assertNotIn(('L0', 'R9'), cands)
        self.assertEqual(int((annotated.positive == False).sum()), 1)  # noqa: E712


class Metrics(unittest.TestCase):
    def test_auroc_aupr(self):
        self.assertEqual(le.auroc([.9, .8], [.1, .8]), .875)
        self.assertAlmostEqual(le.aupr([.9, .8], [.1, .7]), 1.0)
        self.assertAlmostEqual(le.aupr([.5], [.5]), .5)  # one tie group

    def test_ranking_and_consistency(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            pos = pd.DataFrame(dict(ligand=['A', 'B'], receptor=['X', 'Y'], ligand_accessions=['a', 'b'],
                                    receptor_accessions=['x', 'y']))
            pos.to_csv(folder/'positives.csv', index=False)
            pd.DataFrame(dict(ligand=['A', 'B', 'A'], receptor=['X', 'Y', 'Z'], positive=[True, True, False],
                              sources=['t', 't', 't'])).to_csv(folder/'pool_annotated.csv', index=False)
            pos2 = pd.concat([pos, pd.DataFrame(dict(ligand=['C'], receptor=['Z'], ligand_accessions=['c'],
                                                     receptor_accessions=['z']))])
            pos2.to_csv(folder/'positives.csv', index=False)
            p = {('a', 'x'): .2, ('a', 'y'): .9, ('a', 'z'): .99, ('b', 'x'): .1, ('b', 'y'): .5, ('b', 'z'): .5,
                 ('c', 'x'): .3, ('c', 'y'): .3, ('c', 'z'): .1}
            rows = [dict(pair_id=f'{q}__{t}', score=v) for (q, t), v in p.items()]
            rows += [dict(pair_id=f'{t}__{q}', score=v) for (q, t), v in p.items()]
            path = folder/'scores.csv'
            pd.DataFrame(rows).to_csv(path, index=False)
            ranks = le.ranking(folder, le.Scores(path)).set_index('ligand')
            self.assertEqual(ranks.loc['A', 'best_rank'], 2)   # Z is annotated elsewhere: removed; Y beats X
            self.assertEqual(ranks.loc['B', 'best_rank'], 2)   # tie with Z counts against
            self.assertEqual(ranks.loc['C', 'best_rank'], 3)
            off = pd.DataFrame(rows)
            off['score'] += 1e-6
            off['truncated'] = False
            off.loc[0, 'truncated'] = True
            off.loc[0, 'score'] = 0.0
            off.to_csv(folder/'official.csv', index=False)
            c = le.consistency(le.Scores(folder/'official.csv'), le.Scores(path), 'toy', 'm')
            self.assertEqual(c['compared'], len(rows) - 1)
            self.assertTrue(c['within_tolerance'])
            self.assertAlmostEqual(c['max_abs_diff_order'], 1e-6)
        self.assertEqual(le.reverse_id('A__B'), 'B__A')
        self.assertEqual(le.reverse_id('IL6__IL6R__LR'), 'IL6__IL6R__RL')


if __name__ == '__main__':
    unittest.main()
