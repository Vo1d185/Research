import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from diagnostics.human_pairs import build
from diagnostics.summarize import load, summarize


def entry(acc, gene, seq, taxon=9606):
    return {'primaryAccession': acc, 'uniProtkbId': f'{gene}_HUMAN', 'entryType': 'UniProtKB reviewed (Swiss-Prot)',
            'genes': [{'geneName': {'value': gene}}], 'organism': {'taxonId': taxon},
            'proteinDescription': {'recommendedName': {'fullName': {'value': gene}}},
            'entryAudit': {'entryVersion': 1, 'sequenceVersion': 1}, 'sequence': {'value': seq}}


ENTRIES = {'LIG': [entry('P1', 'LIG', 'A' * 50)], 'SHORTR': [entry('P2', 'SHORTR', 'C' * 300)],
           'LONGR': [entry('P3', 'LONGR', 'D' * 900)],
           'AMB': [entry('P4', 'AMB', 'E' * 60), entry('P5', 'AMB', 'F' * 70)]}


def fake_fetch(url):
    headers = {'x-uniprot-release': 'test'}
    if '/search?' in url:
        gene = next(g for g in ENTRIES if f'gene_exact%3A{g}%29' in url)
        return {'results': ENTRIES[gene]}, headers
    acc = url.rsplit('/', 1)[1].split('.')[0]
    return next(e for es in ENTRIES.values() for e in es if e['primaryAccession'] == acc), headers


class DiagnosticsTests(unittest.TestCase):
    def test_groups_orders_and_cross(self):
        pairs = pd.DataFrame(dict(ligand=['LIG', 'LIG'], receptor=['SHORTR', 'LONGR']))
        _, jobs, _ = build(pairs, fetch=fake_fetch)
        cog = jobs[jobs.set == 'cognate']
        self.assertEqual(len(cog), 4)  # both orders
        self.assertEqual(dict(zip(cog.receptor, cog.receptor_group)),
                         {'SHORTR': 'receptor_le_800', 'LONGR': 'receptor_gt_800'})
        self.assertFalse((jobs.set == 'cross').any())  # every ligand-receptor combo is cognate
        rl = jobs[(jobs.receptor == 'SHORTR') & (jobs.order == 'RL')].iloc[0]
        self.assertEqual((rl.query_symbol, rl.text_symbol), ('SHORTR', 'LIG'))

    def test_ambiguous_refused_override_accepted(self):
        pairs = pd.DataFrame(dict(ligand=['LIG'], receptor=['AMB']))
        with self.assertRaisesRegex(ValueError, 'AMB: 2 reviewed'):
            build(pairs, fetch=fake_fetch)
        pmap, _, _ = build(pairs, overrides={'AMB': 'P5'}, fetch=fake_fetch)
        self.assertEqual(pmap.set_index('gene').loc['AMB', 'accession'], 'P5')
        with self.assertRaisesRegex(ValueError, 'not in the pair list'):
            build(pairs, overrides={'XYZ': 'P5'}, fetch=fake_fetch)

    def test_summary_uses_mean_of_orders_and_excludes_known_cross(self):
        ENTRIES['LIG2'] = [entry('P6', 'LIG2', 'G' * 40)]
        pairs = pd.DataFrame(dict(ligand=['LIG', 'LIG2'], receptor=['SHORTR', 'LONGR']))
        _, jobs, _ = build(pairs, known={('LIG2', 'SHORTR'): 'documented'}, fetch=fake_fetch)
        with tempfile.TemporaryDirectory() as d:
            run = Path(d)
            jobs.drop(columns=['query', 'text']).to_csv(run/'jobs.csv', index=False)
            jobs[['pair_id', 'query', 'text']].to_csv(run/'plm_input.csv', index=False)
            score = {'LIG__SHORTR__LR': .9, 'LIG__SHORTR__RL': .7, 'LIG2__LONGR__LR': .2, 'LIG2__LONGR__RL': .4,
                     'LIG__LONGR__LR': .1, 'LIG__LONGR__RL': .1, 'LIG2__SHORTR__LR': .99, 'LIG2__SHORTR__RL': .99}
            pd.DataFrame(dict(pair_id=list(score), score=list(score.values()))).to_csv(run/'scores_x.csv', index=False)
            from proteins import digest
            (run/'scores_x.csv.meta.json').write_text(json.dumps(dict(status='complete',
                input_sha256=digest(run/'plm_input.csv'), output_sha256=digest(run/'scores_x.csv'),
                checkpoint_name='org/ckpt', checkpoint_revision='r')))
            s = load(run)
            self.assertAlmostEqual(s.set_index(['ligand', 'receptor']).loc[('LIG', 'SHORTR'), 'score'], .8)
            groups, ranks = summarize(s)
            le = groups.set_index('receptor_group').loc['receptor_le_800']
            self.assertEqual(le.n_cross, 0)  # the only short-receptor cross pair is documented, so excluded
            self.assertEqual(ranks.set_index('ligand').loc['LIG2', 'rank_among_receptors'], 2)


if __name__ == '__main__':
    unittest.main()
