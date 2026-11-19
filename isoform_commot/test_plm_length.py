"""plm_length.py (reports/plm_length_01): data preparation, padding, merge and the pre-registered decisions."""
import argparse
import json
from pathlib import Path
import tempfile
import unittest
import zlib

import numpy as np
import pandas as pd

import plm_length as pl

AA = 'ACDEFGHIKLMNPQRSTVWY'


def toy(folder, n_prot=60, n_pairs=400, seed=0):
    """Tiny Bernett-like test file and V12-like training file (sequences, labels)."""
    rng = np.random.default_rng(seed)
    lengths = np.r_[rng.integers(60, 400, n_prot - 15), rng.integers(900, 1400, 15)]
    prots = {f'P{i:03d}': ''.join(rng.choice(list(AA), size=n)) for i, n in enumerate(lengths)}
    names = list(prots)
    rows, seen = [], set()
    while len(rows) < n_pairs:
        a, b = rng.choice(names, 2, replace=False)
        if frozenset((a, b)) in seen:
            continue
        seen.add(frozenset((a, b)))
        rows.append(dict(Uniprot_a=a, Uniprot_b=b, label=int(len(rows) % 2 == 0), query=prots[a], text=prots[b]))
    test = pd.DataFrame(rows)
    test.to_csv(folder/'test.csv', index=False)
    train_rows = [dict(query=test['query'][0], text=test['text'][0], label=1),       # pair 0 seen in training
                  dict(query=test['text'][1], text=test['query'][1], label=1),       # pair 1 seen (other order), train pos / test neg
                  dict(query=prots[names[0]], text=prots[names[1]], label=0)]
    train_rows += [dict(query=prots[n], text=prots[names[(i + 3) % 20]], label=0) for i, n in enumerate(names[:20])]
    root = folder/'training'; (root/'STRING_V12_TrainingSet').mkdir(parents=True)
    pd.DataFrame(train_rows).to_csv(root/pl.TRAIN_REL, index=False)
    return test, prots, root


def prep(folder, **kw):
    test, prots, root = toy(folder)
    args = argparse.Namespace(test_file=str(folder/'test.csv'), training_root=str(root), output_dir=str(folder/'runs'/'data'),
                              limit=0, e2_n=kw.get('e2_n', 20), chunk=kw.get('chunk', 300), skip_sha=True)
    pl.prepare(args)
    return test, prots, folder/'runs'


class Prepare(unittest.TestCase):
    def test_overlap_padding_jobs(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            test, prots, runs = prep(d)
            pairs = pd.read_csv(runs/'data'/'pairs.csv')
            self.assertEqual(list(pairs.overlap[:2]), ['pair_seen', 'pair_seen'])
            train = pd.read_csv(d/'training'/pl.TRAIN_REL)
            seen = set(map(frozenset, zip(train['query'], train['text'])))
            expect = [frozenset((q, t)) in seen for q, t in zip(test['query'], test['text'])]
            self.assertEqual(list(~pairs.primary), expect)
            info = json.loads((runs/'data'/'dataset.json').read_text())
            self.assertEqual(info['excluded_pair_seen'], sum(expect))
            pos_seen = set(map(frozenset, zip(train['query'][train.label == 1], train['text'][train.label == 1])))
            self.assertEqual(info['excluded_train_pos_test_neg'],
                             sum(frozenset((q, t)) in pos_seen and l == 0 for q, t, l in zip(test['query'], test['text'], test.label)))
            self.assertEqual(int((pairs.e2 & (pairs.label == 1)).sum()), 20)
            self.assertTrue((pairs.tokens[pairs.e2] <= pl.E2_MAX_TOKENS).all())
            jobs = pd.concat([pd.read_csv(f) for f in sorted((runs/'data'/'jobs').glob('chunk_*.csv'))])
            self.assertTrue(jobs.pair_id.is_unique)
            self.assertEqual(len(jobs), 2 * int(pairs.primary.sum()) + 2 * 40 * len(pl.TARGETS) * len(pl.PADS))
            pos = set(map(frozenset, zip(test.Uniprot_a[test.label == 1], test.Uniprot_b[test.label == 1])))
            for jid, q, t in jobs[jobs.pair_id.str.contains('\\|(?:shuffle|foreign)')].itertuples(index=False):
                uid, variant, order = jid.split('|')
                row = pairs[pairs.uid == uid].iloc[0]
                a, b = prots[row.Uniprot_a], prots[row.Uniprot_b]
                target = int(''.join(c for c in variant if c.isdigit()))
                self.assertEqual(len(q) + len(t) + 3, target)
                qa, qb = (q, t) if order == 'AB' else (t, q)
                self.assertTrue(qa.startswith(a) and qb.startswith(b))
                if variant.startswith('foreign'):     # padding never comes from a partner: check a 30-aa window
                    window = qb[len(b):len(b) + 30]
                    partners = {x for p in pos for x in p if row.Uniprot_a in p or row.Uniprot_b in p}
                    self.assertFalse(any(window in prots[x] for x in partners | {row.Uniprot_a, row.Uniprot_b}
                                         if len(window) == 30))
            again = d/'again'; again.mkdir()
            pl.prepare(argparse.Namespace(test_file=str(d/'test.csv'), training_root=str(d/'training'),
                                          output_dir=str(again/'data'), limit=0, e2_n=20, chunk=300, skip_sha=True))
            for f in sorted((runs/'data'/'jobs').glob('chunk_*.csv')):
                self.assertEqual(f.read_bytes(), (again/'data'/'jobs'/f.name).read_bytes())   # deterministic


def fake_scores(runs, signal):
    """Write plm_score.py-like outputs: score = sigmoid(label-dependent signal(tokens) + noise)."""
    pairs = pd.read_csv(runs/'data'/'pairs.csv').set_index('uid')
    for chunk in sorted((runs/'data'/'jobs').glob('chunk_*.csv')):
        j = pd.read_csv(chunk)
        uid = j.pair_id.str.split('|').str[0]
        tok = j['query'].str.len() + j['text'].str.len() + 3
        lab = pairs.loc[uid, 'label'].to_numpy()
        rng = np.random.default_rng(zlib.crc32(chunk.name.encode()))
        z = (2 * lab - 1) * signal(tok.to_numpy()) + rng.normal(size=len(j))
        out = runs/'scores'/chunk.name; out.parent.mkdir(exist_ok=True, parents=True)
        pd.DataFrame(dict(pair_id=j.pair_id, score=1 / (1 + np.exp(-z)), n_tokens=tok)).to_csv(out, index=False)
        pl.write_json(Path(str(out) + '.meta.json'), dict(status='complete', input_sha256=pl.digest(chunk),
                                                          checkpoint_sha256='x'))


class MergeEvaluate(unittest.TestCase):
    def run_case(self, signal):
        d = Path(tempfile.mkdtemp())
        _, _, runs = prep(d, e2_n=20)
        fake_scores(runs, signal)
        pl.merge(argparse.Namespace(runs=str(runs), output=str(runs/'scores.csv')))
        pl.evaluate(argparse.Namespace(runs=str(runs), scores=str(runs/'scores.csv'), output_dir=str(d/'res'),
                                       lr_scores=None, n_boot=200))
        return json.loads((d/'res'/'verdict.json').read_text()), runs

    def test_merge_averages_orders(self):
        _, runs = self.run_case(lambda t: np.full(len(t), 1.0))
        s = pd.read_csv(runs/'scores.csv')
        self.assertTrue(np.allclose(s.score, (s.score_ab + s.score_ba) / 2))

    def test_merge_detects_truncation(self):
        d = Path(tempfile.mkdtemp())
        _, _, runs = prep(d)
        fake_scores(runs, lambda t: np.ones(len(t)))
        f = sorted((runs/'scores').glob('chunk_*.csv'))[0]
        s = pd.read_csv(f); s.loc[0, 'n_tokens'] -= 1; s.to_csv(f, index=False)
        with self.assertRaises(ValueError):
            pl.merge(argparse.Namespace(runs=str(runs), output=str(runs/'scores.csv')))

    def test_length_effect_detected(self):
        v, _ = self.run_case(lambda t: np.where(t <= 1603, 3.0, 0.3))
        self.assertIn(v['E1'], ('long_clearly_worse', 'long_slightly_worse'))
        self.assertEqual(v['E2'], 'length_itself_lowers_discrimination')


class Decisions(unittest.TestCase):
    def test_auroc(self):
        self.assertEqual(pl.auroc([3, 4], [1, 2]), 1.0); self.assertEqual(pl.auroc([1], [1]), 0.5)

    def test_e1_rules(self):
        st = {k: dict(delta=0.1) for k in pl.STRATA}
        self.assertEqual(pl.e1_decision(dict(delta=0.1, ci95=[0.05, 0.15]), st), 'long_clearly_worse')
        st2 = dict(st, none_seen=dict(delta=-0.01))
        self.assertEqual(pl.e1_decision(dict(delta=0.1, ci95=[0.05, 0.15]), st2), 'long_slightly_worse')
        self.assertEqual(pl.e1_decision(dict(delta=0.01, ci95=[-0.02, 0.04]), st), 'no_material_difference')
        self.assertEqual(pl.e1_decision(dict(delta=0.01, ci95=[-0.08, 0.1]), st), 'uncertain')

    def test_e2_rules(self):
        t = pd.DataFrame(dict(target=[2400, 2400], pad=list(pl.PADS), delta=[0.1, 0.08], ci95_low=[0.05, 0.02], ci95_high=[0.15, 0.12]))
        self.assertEqual(pl.e2_decision(t), 'length_itself_lowers_discrimination')
        t2 = t.assign(delta=[0.0, 0.01], ci95_low=[-0.02, -0.01], ci95_high=[0.02, 0.03])
        self.assertEqual(pl.e2_decision(t2), 'length_itself_small_effect')
        self.assertEqual(pl.e2_decision(t.assign(ci95_low=[0.05, -0.01])), 'uncertain')

    def test_lr_replication(self):
        rng = np.random.default_rng(3)
        n = 400
        rows = []
        for i in range(n):
            pos = i % 5 == 0
            long_ = i % 2 == 0
            L = 1200 if long_ else 300
            signal = (0.4 if long_ else 2.5) * (1 if pos else -1)
            rows.append(dict(negative='positive' if pos else 'random', training=['proteins_seen', 'protein_unseen', 'seen_positive'][i % 3],
                             ligand_lengths=str(L), receptor_lengths=str(L // 2), score=1 / (1 + np.exp(-(signal + rng.normal())))))
        path = Path(tempfile.mkdtemp())/'ha.csv'
        pd.DataFrame(rows).to_csv(path, index=False)
        res = pl.lr_replication(path, np.random.default_rng(1))
        self.assertGreater(res['overall']['delta'], 0.1)
        kept = sum(r['training'] != 'seen_positive' for r in rows)
        self.assertEqual(sum(res['overall']['n'].values()), kept)

    def test_allocate(self):
        self.assertEqual(pl.allocate({'a': 1, 'b': 1, 'c': 1}, 2), {'a': 1, 'b': 1, 'c': 0})
        self.assertEqual(sum(pl.allocate({'a': 5, 'b': 3, 'c': 2}, 7).values()), 7)


if __name__ == '__main__':
    unittest.main()
