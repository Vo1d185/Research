"""Toy tests for stage2_ppi.py: subunit expansion, homophilic pairs, job de-duplication, import rules, designs."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

import stage2_ppi as sp

SEQ = {'L1': 'ACDEFGHIKL', 'R1': 'MNPQRSTVWY' * 2, 'R2': 'ACDACDACDACD', 'L2': 'KLMNKLMN',
       'G1': 'WYWYWYWYWY', 'G2': 'PQPQPQPQ', 'L3': 'EEEEKKKKEE'}
ACC = {g: f'P{i:05d}' for i, g in enumerate(sorted(SEQ))}


def toy(tmp, unresolved=('L4',)):
    inp = Path(tmp)/'input'; inp.mkdir()
    pd.DataFrame([('L1', 'R1', 'diffusion'), ('L1', 'R2', 'diffusion'), ('L2', 'CX', 'contact'),
                  ('L3', 'L3', 'contact'), ('L4', 'R1', 'diffusion')],
                 columns=['ligand', 'receptor', 'type']).to_csv(inp/'lr.csv', index=False)
    ents = {'L1': 'L1', 'R1': 'R1', 'R2': 'R2', 'L2': 'L2', 'CX': 'G1;G2', 'L3': 'L3', 'L4': 'L4'}
    pd.DataFrame(dict(entity=list(ents), subunits=list(ents.values()))).to_csv(inp/'entities.csv', index=False)
    prot = Path(tmp)/'proteins'; prot.mkdir()
    pd.DataFrame(dict(lr_symbol=list(ACC), data_symbol=list(ACC), accession=list(ACC.values()))).to_csv(
        prot/'protein_map.csv', index=False)
    (prot/'proteins.fasta').write_text(''.join(f'>{ACC[g]}|{g}|{g}\n{SEQ[g]}\n' for g in ACC))
    (prot/'proteins.json').write_text(json.dumps(dict(status='complete', unresolved=list(unresolved),
        fasta_sha256=hashlib.sha256((prot/'proteins.fasta').read_bytes()).hexdigest())))
    return inp, prot


def fake_scores(tmp, jobs_dir, fn, sha=sp.EXPECTED_CHECKPOINT_SHA256):
    jb = pd.read_csv(Path(jobs_dir)/'jobs.csv')
    sc = pd.DataFrame(dict(pair_id=jb.pair_id, score=[fn(r) for r in jb.itertuples()], n_tokens=jb.n_tokens))
    path = Path(tmp)/'scores.csv'; sc.to_csv(path, index=False)
    Path(str(path) + '.meta.json').write_text(json.dumps(dict(
        status='complete', input_sha256=sp.digest(Path(jobs_dir)/'plm_input.csv'), checkpoint_sha256=sha)))
    return path


class Cli(unittest.TestCase):
    def test_repeated_input_flags_accumulate(self):
        a = sp.parser().parse_args(['jobs', '--input', 'A', '--input', 'B', 'C', '--proteins', 'P', '--output', 'O'])
        self.assertEqual(a.input, ['A', 'B', 'C'])


class Panel(unittest.TestCase):
    def test_genes_and_jobs(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp, prot = toy(tmp)
            sp.genes(argparse.Namespace(input=[inp], output=Path(tmp)/'genes.csv'))
            g = pd.read_csv(Path(tmp)/'genes.csv')
            self.assertEqual(set(g.symbol), {'L1', 'R1', 'R2', 'L2', 'G1', 'G2', 'L3', 'L4'})
            self.assertTrue(bool(g.set_index('symbol').loc['L3', 'as_ligand'] and g.set_index('symbol').loc['L3', 'as_receptor']))
            sp.jobs(argparse.Namespace(input=[inp], proteins=prot, output=Path(tmp)/'jobs', limit_pairs=0,
                                       include_longest=False, seed=0))
            ent = pd.read_csv(Path(tmp)/'jobs'/'entity_pairs.csv').set_index(['ligand', 'receptor'])
            self.assertFalse(ent.loc[('L4', 'R1'), 'complete'])                  # unmapped subunit -> no score
            self.assertTrue(ent.loc[('L3', 'L3'), 'homophilic'])
            self.assertEqual(len(ent.loc[('L2', 'CX'), 'subunit_pair_keys'].split(';')), 2)
            jb = pd.read_csv(Path(tmp)/'jobs'/'jobs.csv')
            # unique subunit pairs: L1-R1, L1-R2, L2-G1, L2-G2, L3-L3  -> 5 keys x 2 orders
            self.assertEqual(jb.key.nunique(), 5); self.assertEqual(len(jb), 10)
            self.assertTrue((jb.n_tokens == jb['query'].str.len() + jb['text'].str.len() + 3).all())


class Import(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); t = self.tmp.name
        self.inp, prot = toy(t)
        sp.jobs(argparse.Namespace(input=[self.inp], proteins=prot, output=Path(t)/'jobs', limit_pairs=0,
                                   include_longest=False, seed=0))
        self.jobs = Path(t)/'jobs'

    def tearDown(self):
        self.tmp.cleanup()

    def test_rules(self):
        t = self.tmp.name
        g1, g2 = ACC['G1'], ACC['G2']
        def fn(r):          # complex: G1 subunit 0.2 (both orders), G2 subunit 0.6 / 0.8 -> mean 0.7, max 0.7
            if g2 in r.key:
                return 0.6 if r.order == 'AB' else 0.8
            return 0.2 if g1 in r.key else 0.5
        path = fake_scores(t, self.jobs, fn)
        sp.import_scores(argparse.Namespace(jobs=self.jobs, scores=path, output=Path(t)/'ppi', exploratory=False))
        ppi = pd.read_csv(Path(t)/'ppi'/'ppi_scores.csv').set_index(['ligand', 'receptor'])
        self.assertAlmostEqual(ppi.loc[('L2', 'CX'), 'score'], 0.7)
        self.assertTrue(np.isnan(ppi.loc[('L4', 'R1'), 'score']))
        self.assertAlmostEqual(ppi.loc[('L1', 'R1'), 'score'], 0.5)

    def test_refusals(self):
        t = self.tmp.name
        path = fake_scores(t, self.jobs, lambda r: 0.5, sha='other')
        with self.assertRaises(ValueError):
            sp.import_scores(argparse.Namespace(jobs=self.jobs, scores=path, output=Path(t)/'x', exploratory=False))
        sc = pd.read_csv(path).iloc[:-1]; sc.to_csv(path, index=False)
        with self.assertRaises(ValueError):
            sp.import_scores(argparse.Namespace(jobs=self.jobs, scores=path, output=Path(t)/'y', exploratory=True))


class Designs(unittest.TestCase):
    def test_normalisation(self):
        lr = pd.DataFrame([('A', 'X', 'diffusion'), ('A', 'Y', 'diffusion'), ('A', 'Z', 'diffusion'),
                           ('B', 'X', 'diffusion'), ('C', 'Y', 'contact'), ('C', 'Z', 'contact')],
                          columns=['ligand', 'receptor', 'type'])
        ppi = lr.assign(score=[0.9, 0.09, 0.01, 0.5, np.nan, 0.3])
        d = sp.design_table(lr, ppi).set_index(['ligand', 'receptor'])
        self.assertAlmostEqual(d.loc[('A', 'X'), 's_norm_ligand'], 0.9)
        self.assertAlmostEqual(d.loc[pd.IndexSlice['A', :], 's_norm_ligand'].sum(), 1.0)
        self.assertEqual(d.loc[('B', 'X'), 's_norm_ligand'], 1.0)                  # single receptor -> 1
        self.assertTrue(d.loc[pd.IndexSlice['C', :], 's_norm_ligand'].isna().all())   # missing score -> group NaN
        self.assertAlmostEqual(d.loc[('A', 'Y'), 's_uniform_ligand'], 1 / 3)
        self.assertAlmostEqual(d.loc[('A', 'X'), 's_norm_receptor'], 0.9 / 1.4)
        s = sp.describe(d.reset_index())
        self.assertEqual(s['ligands_k1'], 1); self.assertEqual(s['ligands_with_missing_score'], 1)


class Cross(unittest.TestCase):
    def test_sampling_rules(self):
        lr = pd.DataFrame([('L1', 'R1', 'diffusion'), ('L1', 'R2', 'diffusion'), ('L4', 'R1', 'diffusion'),
                           ('L2', 'CX', 'contact'), ('L3', 'L3', 'contact'), ('L5', 'R9', 'contact')],
                          columns=['ligand', 'receptor', 'type'])
        sub = {'L1': ['L1'], 'R1': ['R1'], 'R2': ['R2'], 'L4': ['L4'], 'L2': ['L2'], 'CX': ['G1', 'L5'],
               'L3': ['L3'], 'L5': ['L5'], 'R9': ['R9']}
        cand_all, n = sp.cross_pairs(lr, sub, set(), 0, 0)
        # diffusion: L4-R2; contact: L2-L3, L2-R9, L3-CX, L3-R9, L5-L3 (L5-CX shares subunit L5; L3-L3 self)
        _, n_all = sp.cross_pairs(lr, sub, set(), 1, 0)
        self.assertEqual(n_all, 6)
        _, n_ex = sp.cross_pairs(lr, sub, {('L4', 'R2')}, 1, 0)
        self.assertEqual(n_ex, 5)
        got, _ = sp.cross_pairs(lr, sub, set(), 6, 3)
        self.assertNotIn(('L5', 'CX'), set(zip(got.ligand, got.receptor)))
        self.assertTrue((got.type == got.ligand.map(lr.groupby('ligand').type.first())).all())
        with self.assertRaises(ValueError):
            sp.cross_pairs(lr, sub, set(), 7, 0)

    def test_jobs_and_import_with_cross(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp, prot = toy(tmp)
            ex = Path(tmp)/'ex.csv'; pd.DataFrame(dict(ligand=['L4'], receptor=['R2'])).to_csv(ex, index=False)
            sp.jobs(argparse.Namespace(input=[inp], proteins=prot, output=Path(tmp)/'jobs', limit_pairs=0,
                                       include_longest=False, seed=0, cross_ratio=0.5, cross_seed=1,
                                       cellchat=None, exclude_pairs=str(ex)))
            ent = pd.read_csv(Path(tmp)/'jobs'/'entity_pairs.csv')
            self.assertEqual(int((ent.set == 'cross').sum()), 2)     # 0.5 x 4 complete panel pairs
            path = fake_scores(tmp, Path(tmp)/'jobs', lambda r: 0.9 if r.key.count('P') == 2 else 0.5)
            sp.import_scores(argparse.Namespace(jobs=Path(tmp)/'jobs', scores=path, output=Path(tmp)/'ppi',
                                                exploratory=False))
            rep = json.loads((Path(tmp)/'ppi'/'score_report.json').read_text())
            self.assertIsNotNone(rep['panel_vs_unannotated'])
            self.assertEqual(len(pd.read_csv(Path(tmp)/'ppi'/'ppi_scores.csv')), 5)


if __name__ == '__main__':
    unittest.main()
