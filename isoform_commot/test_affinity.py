"""Unit tests for affinity.py (reports/affinity_01/plan.md). Run: python -m unittest test_affinity"""
import argparse
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

import affinity as af

VERSION = '"# GtoPdb Version: test - published: 2026-01-01"'
L = dict(A='AAAGKWVLTY' * 2, B='CCDEFGHIKL' * 2, C='MNPQRSTVWY' * 2, D='KKWWYYFFLL' * 2)
R = dict(X='GGSAWKTLVY' * 6, Y='PPQMNRSTVC' * 6, Z='LLKKDDEEFF' * 6)


def toy(tmp, extra_rows=(), training_rows=()):
    """GtoPdb-shaped inputs: receptor X with 3 ligands, Y with 3, plus rows the filters must drop."""
    tmp = Path(tmp)
    gd = tmp/'in'; gd.mkdir(parents=True)
    rows, pep = [], []
    def add(target, acc, lig_id, lig, value, units='pKd', endo='true', ltype='Peptide', n=1):
        for _ in range(n):
            rows.append(dict(Target=target, **{'Target UniProt ID': acc, 'Target Species': 'Human',
                        'Ligand ID': lig_id, 'Ligand': 'lig' + str(lig_id), 'Ligand Type': ltype,
                        'Endogenous': endo, 'Affinity Units': units, 'Affinity Median': value,
                        'Affinity High': np.nan, 'Affinity Low': np.nan}))
        if not any(p['Ligand id'] == lig_id for p in pep):
            pep.append({'Ligand id': lig_id, 'Single letter amino acid sequence': lig})
    for i, (lig, v) in enumerate(zip('ABC', [9.0, 8.0, 6.0])):
        add('X', 'P00001', 10 + i, L[lig], v)
    for i, (lig, v) in enumerate(zip('ABD', [7.0, 9.5, 8.5])):
        add('Y', 'P00002', 10 + i if lig != 'D' else 20, L[lig], v, units='pKi')
    add('Z', 'P00003', 30, L['A'], 8.0, endo='false')                 # dropped: not endogenous
    add('Z', 'P00003', 31, L['B'], 8.0, units='pEC50')                # dropped: functional unit
    add('Z', 'P00003', 32, L['C'], 8.0, ltype='Metabolite')           # dropped: not a peptide
    add('W', 'P00001|P00002', 33, L['D'], 8.0)                        # dropped: multi-gene target
    add('Z', 'P00003', 34, 'AAA', 8.0)                                # dropped: sequence too short
    add('Z', 'P00003', 35, 'XAAGKWVLTY', 7.0)                         # kept: pyroglutamate -> Q
    for r in extra_rows:
        add(**r)
    df = pd.DataFrame(rows)
    with (gd/'gtopdb_interactions.csv').open('w') as f:
        f.write(VERSION + '\n'); df.to_csv(f, index=False)
    with (gd/'gtopdb_peptides.csv').open('w') as f:
        f.write(VERSION + '\n'); pd.DataFrame(pep).to_csv(f, index=False)
    up = pd.DataFrame(dict(accession=['P00001', 'P00002', 'P00003'], sequence=[R['X'], R['Y'], R['Z']]))
    up['length'] = up.sequence.str.len()
    up.to_csv(tmp/'uniprot.tsv', sep='\t', index=False)
    tr = pd.DataFrame([{'Receptor Sequence': a, 'Ligand Sequence': b} for a, b in
                       [(R['Z'], L['D'])] + list(training_rows)])
    tr.to_csv(tmp/'train.csv', index=False)
    return gd, tmp/'uniprot.tsv', tmp/'train.csv'


def prepared(tmp, **kw):
    gd, up, tr = toy(tmp, **kw)
    out = Path(tmp)/'data'
    af.prepare(argparse.Namespace(gtopdb_dir=str(gd), uniprot=str(up), training_csv=str(tr),
                                  output_dir=str(out), limit=0))
    return out


class Prepare(unittest.TestCase):
    def test_filters_and_pyroglutamate(self):
        tmp = tempfile.mkdtemp()
        out = prepared(tmp)
        p = pd.read_csv(out/'pairs.csv')
        counts = json.loads((out/'dataset.json').read_text())['counts']
        self.assertEqual(counts['endogenous'], 11)                 # the non-endogenous row is gone
        self.assertEqual(counts['peptide'], 10)
        self.assertEqual(counts['binding_units'], 9)
        self.assertEqual(counts['single_gene_target'], 8)
        self.assertEqual(counts['pyroglutamate_renamed'], 1)
        self.assertEqual(set(p.rec_acc), {'P00001', 'P00002', 'P00003'})
        self.assertTrue(p.lig_seq.map(lambda s: set(s) <= af.AA).all())
        self.assertEqual(len(p), 7)
        self.assertTrue((p.pair_id.str.startswith('human__')).all() and p.pair_id.is_unique)
        q = pd.read_csv(out/'plm_input.csv')                        # ligand first, as in plan section 1
        self.assertEqual(list(q.columns), ['pair_id', 'query', 'text'])
        self.assertEqual(q.iloc[0]['query'], p.iloc[0].lig_seq)

    def test_median_of_repeats(self):
        tmp = tempfile.mkdtemp()
        rows = [dict(target='X', acc='P00001', lig_id=40, lig=L['D'], value=v) for v in (6.0, 8.0, 7.0)]
        out = prepared(tmp, extra_rows=rows)
        p = pd.read_csv(out/'pairs.csv')
        row = p[p.pair_id == 'human__P00001__40'].iloc[0]
        self.assertEqual((row.value, row.n_records), (7.0, 3))

    def test_leakage_excludes_only_both_sides(self):
        tmp = tempfile.mkdtemp()
        out = prepared(tmp, training_rows=[(R['X'], L['A'])])       # receptor X with ligand A was trained on
        p, ex = pd.read_csv(out/'pairs.csv'), pd.read_csv(out/'excluded.csv')
        self.assertEqual(list(ex.pair_id), ['human__P00001__10'])
        kept = p[p.rec_acc == 'P00001']
        self.assertTrue(kept.rec_seen.all())                        # receptor seen, pair not seen -> kept
        self.assertEqual(int(p[p.pair_id == 'human__P00002__10'].lig_seen.iloc[0]), 1)


class Metrics(unittest.TestCase):
    def table(self, bind, v12=None, length=None):
        n = len(bind)
        return pd.DataFrame(dict(species='Human', rec_acc='P1', lig_id=range(n), value=np.arange(n, dtype=float),
                                 bindpred=bind, v12=v12 if v12 is not None else np.zeros(n),
                                 length=length if length is not None else np.arange(n, dtype=float),
                                 units='pKd', lig_len=10, lig_seen=False, rec_seen=False))

    def test_group_rho_and_sign(self):
        g = af.group_rho(self.table([0, 1, 2, 3]), ['species', 'rec_acc'], 'lig_id')
        self.assertAlmostEqual(g.bindpred.iloc[0], 1.0)
        g2 = af.group_rho(self.table([3, 2, 1, 0], length=[3, 2, 1, 0]), ['species', 'rec_acc'], 'lig_id')
        self.assertAlmostEqual(g2.bindpred.iloc[0], -1.0)
        self.assertAlmostEqual(af.rho_bars(g2)['length'], 1.0)      # baseline takes the favourable sign
        self.assertAlmostEqual(af.rho_bars(g2)['bindpred'], -1.0)
        small = af.group_rho(self.table([0, 1]), ['species', 'rec_acc'], 'lig_id')
        self.assertEqual(len(small), 0)                             # fewer than MIN_GROUP members

    def test_decide_thresholds(self):
        gs = pd.DataFrame(dict(group=[f'g{i}' for i in range(12)], n=3, value_range=2.0, informative_pairs=3,
                               bindpred=[0.9] * 12, v12=[0.1] * 12, length=[0.05] * 12))
        d = af.decide(gs)
        self.assertEqual((d['A1_above_zero'], d['A2_effect_size'], d['A3_above_baselines']), (True, True, True))
        self.assertEqual(d['verdict'], '可作为 S 进入 OT')
        weak = gs.assign(bindpred=[0.1] * 12)
        dw = af.decide(weak)
        self.assertEqual((dw['A1_above_zero'], dw['A2_effect_size']), (True, False))
        self.assertEqual(af.decide(gs.head(5))['undetermined'], True)


try:
    from catboost import CatBoostRegressor
    HAS_CATBOOST = True
except ImportError:                       # predict runs in --predict-env; see reports/affinity_01/README.md
    HAS_CATBOOST = False


@unittest.skipUnless(HAS_CATBOOST, 'catboost is not installed in this environment')
class EndToEnd(unittest.TestCase):
    def test_predict_and_evaluate(self):
        tmp = Path(tempfile.mkdtemp())
        out = prepared(tmp)
        pairs = pd.read_csv(out/'pairs.csv')
        seqs = sorted(set(pairs.lig_seq) | set(pairs.rec_seq))
        rng = np.random.default_rng(0)
        vec = {s: rng.normal(size=4) for s in seqs}
        np.savez_compressed(tmp/'emb.npz', sequences=np.asarray(seqs, object),
                            vectors=np.vstack([vec[s] for s in seqs]), meta='{}')
        X = np.vstack([np.r_[vec[l], vec[r]] for l, r in zip(pairs.lig_seq, pairs.rec_seq)])
        m = CatBoostRegressor(iterations=400, depth=2, learning_rate=0.3, verbose=0)
        m.fit(X, -pairs.value.to_numpy())                            # toy model: learns the labels
        m.save_model(str(tmp/'model.cbm'), format='cbm')
        af.predict(argparse.Namespace(pairs=str(out/'pairs.csv'), embeddings=str(tmp/'emb.npz'),
                                      model=str(tmp/'model.cbm'), output=str(tmp/'bindpred.csv')))
        bp = pd.read_csv(tmp/'bindpred.csv')
        self.assertEqual(len(bp), len(pairs))
        np.testing.assert_allclose(bp.pkd_hat, -bp.log10_kd_pred)
        pd.DataFrame(dict(pair_id=pairs.pair_id, score=rng.random(len(pairs)), n_tokens=10)) \
            .to_csv(tmp/'v12.csv', index=False)
        af.evaluate(argparse.Namespace(pairs=str(out/'pairs.csv'), bindpred=str(tmp/'bindpred.csv'),
                                       v12=str(tmp/'v12.csv'), output=str(tmp/'summary')))
        dec = json.loads((tmp/'summary'/'decisions.json').read_text())
        self.assertTrue((tmp/'summary'/'affinity01_figure.png').exists())
        self.assertIn('verdict', dec['primary_receptor_groups'] | {'verdict': None})
        self.assertEqual(dec['primary_receptor_groups']['undetermined'], True)   # only 2 toy groups
        g = pd.read_csv(tmp/'summary'/'groups_receptor.csv')
        self.assertEqual(len(g), 2)
        self.assertTrue((g.bindpred > 0.5).all())                    # the toy model recovers the ranking
        scored = pd.read_csv(tmp/'summary'/'pairs_scored.csv')
        self.assertNotIn('rec_seq', scored.columns)
        np.testing.assert_allclose(scored.log10_kd_exp, -scored.value)


if __name__ == '__main__':
    unittest.main()
