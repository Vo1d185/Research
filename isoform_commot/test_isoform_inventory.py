import argparse
import json
import random
import tempfile
import unittest
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

import isoform_inventory as inv


def feature(kind, start, end, description=''):
    return dict(type=kind, description=description, location=dict(start=dict(value=start), end=dict(value=end)))


def entry(accession, gene, length, features=(), comments=(), sequence=None, mgi=None):
    return dict(primaryAccession=accession, uniProtkbId=accession + '_MOUSE', entryType='UniProtKB reviewed (Swiss-Prot)',
                genes=[dict(geneName=dict(value=gene))], organism=dict(taxonId=10090),
                proteinDescription=dict(recommendedName=dict(fullName=dict(value=gene + ' protein'))),
                entryAudit=dict(entryVersion=1, sequenceVersion=1), sequence=dict(value=sequence or 'A' * length),
                features=list(features), comments=list(comments),
                uniProtKBCrossReferences=[dict(database='MGI', id=mgi)] if mgi else [])


def protein(n, seed):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(n))


GX = protein(300, 1)            # alpha-like reference entry of gene Gx
GY = protein(400, 2)


def alt_products(n):
    return dict(commentType='ALTERNATIVE PRODUCTS',
                isoforms=[dict(name=dict(value=str(i + 1)), isoformIds=[f'X-{i + 1}'],
                               isoformSequenceStatus='Displayed' if i == 0 else 'Described') for i in range(n)])


ENTRIES = {
    # receptor: variable region inside the extracellular Ig domain -> candidate
    'Rec1': [entry('P1', 'Rec1', 300, [feature('Topological domain', 20, 150, 'Extracellular'),
                                       feature('Transmembrane', 151, 171), feature('Topological domain', 172, 300, 'Cytoplasmic'),
                                       feature('Domain', 40, 120, 'Ig-like'), feature('Alternative sequence', 60, 80, 'in isoform 2')],
                   [alt_products(2)])],
    # receptor: variable region only in the cytoplasmic tail -> not a candidate
    'Rec2': [entry('P2', 'Rec2', 300, [feature('Topological domain', 20, 150, 'Extracellular'),
                                       feature('Topological domain', 172, 300, 'Cytoplasmic'),
                                       feature('Alternative sequence', 250, 260, 'in isoform 2')], [alt_products(2)])],
    # secreted ligand: no topology; variable region after the signal peptide -> extracellular
    'Lig1': [entry('P3', 'Lig1', 200, [feature('Signal', 1, 25), feature('Alternative sequence', 100, 120, 'in isoform 3')],
                   [alt_products(3), dict(commentType='SUBCELLULAR LOCATION',
                                          subcellularLocations=[dict(location=dict(value='Secreted'))])])],
    # two reviewed entries (alpha/beta-like), no alternative sequences
    'Nrxn1': [entry('Q9CS84', 'Nrxn1', 1514), entry('P0DI97', 'Nrxn1', 468)],
    'Sub2': [],
    # Nrxn2-like: two reviewed entries, no alternative sequences; the short entry replaces the reference's
    # residues 1-200 (extracellular, Laminin G-like domain) by 20 own residues -> candidate by the several-entries rule
    'Gx': [entry('G1', 'Gx', 300, [feature('Topological domain', 1, 250, 'Extracellular'),
                                   feature('Domain', 50, 150, 'Laminin G-like 1')], [alt_products(2)], GX, 'MGI:1'),
           entry('G2', 'Gx', 120, [], [alt_products(2)], 'M' + protein(19, 3) + GX[200:], 'MGI:1')],
    # second entry of Gy has another primary gene name, found only through the shared MGI id; it lacks 10
    # cytoplasmic residues -> not a candidate
    'Gy': [entry('Y1', 'Gy', 400, [feature('Topological domain', 1, 200, 'Extracellular'),
                                   feature('Topological domain', 221, 400, 'Cytoplasmic')], [], GY, 'MGI:2')],
}
XREF_ONLY = [entry('Y2', 'Gyb', 390, [], [], GY[:300] + GY[310:], 'MGI:2')]


def fake_fetch(url):
    url = unquote(url)
    if 'gene_exact:' in url:
        gene = url.split('gene_exact:')[1].split(')')[0]
        return dict(results=ENTRIES.get(gene, [])), {'x-uniprot-release': '2026_05'}
    mgi = url.split('xref:"mgi-')[1].split('"')[0]
    pool = [e for es in ENTRIES.values() for e in es] + XREF_ONLY
    return dict(results=[e for e in pool if dict(database='MGI', id=mgi) in e['uniProtKBCrossReferences']]), \
        {'x-uniprot-release': '2026_05'}


class Inventory(unittest.TestCase):
    def test_polite_retries_rate_limit(self):
        from urllib.error import HTTPError
        calls = []

        def flaky(url):
            calls.append(url)
            if len(calls) < 3:
                raise HTTPError(url, 429, 'Too Many Requests', {'Retry-After': '0'}, None)
            return 'ok'
        self.assertEqual(inv.polite(flaky)('u'), 'ok')
        self.assertEqual(len(calls), 3)

    def test_flags_candidates_and_known_example_first(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root/'in').mkdir()
            pd.DataFrame(dict(ligand=['Lig1', 'Nrxn1', 'Lig1', 'Gx'], receptor=['Rec1', 'Rec2', 'CX', 'Gy'],
                              type=['diffusion', 'contact', 'diffusion', 'contact'])).to_csv(root/'in'/'lr.csv',
                                                                                             index=False)
            pd.DataFrame(dict(entity=['Lig1', 'Nrxn1', 'Rec1', 'Rec2', 'CX', 'Gx', 'Gy'],
                              subunits=['Lig1', 'Nrxn1', 'Rec1', 'Rec2', 'Rec1;Sub2', 'Gx', 'Gy'])).to_csv(
                root/'in'/'entities.csv', index=False)
            (root/'in'/'preparation.json').write_text(json.dumps(dict(sample='synthetic')))
            inv.run(argparse.Namespace(input=str(root/'in'), output=str(root/'out'), workers=2, previous=None),
                    fetch=fake_fetch)
            genes = pd.read_csv(root/'out'/'genes.csv').set_index('gene')
            self.assertTrue(genes.loc['Rec1', 'candidate'] and genes.loc['Rec1', 'any_extracellular'])
            self.assertTrue(genes.loc['Rec1', 'any_domain'])
            self.assertFalse(genes.loc['Rec2', 'candidate'])                 # cytoplasmic only
            self.assertTrue(genes.loc['Lig1', 'candidate'] and genes.loc['Lig1', 'any_extracellular'])
            self.assertEqual(genes.loc['Nrxn1', 'reviewed_entries'], 2)
            self.assertEqual(genes.loc['Sub2', 'reviewed_entries'], 0)
            self.assertEqual(set(genes.loc['Rec1', 'lr_pairs'].split(';')), {'Lig1->Rec1', 'Lig1->CX'})  # complex subunit
            candidates = pd.read_csv(root/'out'/'candidates.csv')
            self.assertEqual(candidates.gene.iloc[0], 'Nrxn1')                # known example first
            self.assertEqual(set(candidates.gene), {'Nrxn1', 'Rec1', 'Lig1', 'Gx'})
            self.assertTrue(genes.loc['Gx', 'candidate_several_entries'])
            self.assertFalse(genes.loc['Gx', 'candidate_alternative_sequence'])   # missed by the inventory 01 rule
            self.assertEqual(genes.loc['Gx', 'candidate_rule'], 'several_entries')
            self.assertEqual(genes.loc['Gx', 'reference_entry'], 'G1')
            self.assertTrue(genes.loc['Gx', 'between_extracellular'] and genes.loc['Gx', 'between_domain'])
            self.assertEqual((genes.loc['Gy', 'reviewed_entries'], genes.loc['Gy', 'entries_by_mgi_only']), (2, 1))
            self.assertFalse(genes.loc['Gy', 'candidate'])                     # difference only in the cytoplasm
            self.assertEqual(genes.loc['Gy', 'merged_isoforms'], 2)            # Y1-1, Y2-1
            regions = pd.read_csv(root/'out'/'variable_regions.csv')
            gx = regions[regions.gene == 'Gx']
            self.assertEqual(gx[['source', 'start', 'end', 'kind']].values.tolist(),
                             [['between_entries', 1, 200, 'replaced']])
            meta = json.loads((root/'out'/'inventory.json').read_text())
            self.assertEqual(meta['uniprot_releases'], ['2026_05'])
            self.assertEqual(meta['genes_with_several_entries'], ['Gx', 'Gy', 'Nrxn1'])
            self.assertEqual(meta['entries_found_by_mgi_only'], ['Y2'])
            self.assertEqual(meta['candidates_several_entries_only'], ['Gx'])

    def test_entry_differences(self):
        ref = protein(100, 4)
        regions, shared = inv.entry_differences(ref, ref[:40] + 'WWWWW' + ref[40:90])     # insertion, end missing
        self.assertEqual(shared, 90)
        self.assertEqual([(r['start'], r['end'], r['kind'], r['other_residues']) for r in regions],
                         [(40, 41, 'inserted', 5), (91, 100, 'absent', 0)])
        regions, _ = inv.entry_differences(ref, ref[:30] + ref[35:])                     # short internal deletion
        self.assertEqual([(r['start'], r['end'], r['kind']) for r in regions], [(31, 35, 'absent')])
        regions, shared = inv.entry_differences(ref, ref)
        self.assertEqual((regions, shared), ([], 100))


if __name__ == '__main__':
    unittest.main()
