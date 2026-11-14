import argparse
import json
import tempfile
import unittest
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

import isoform_inventory as inv


def feature(kind, start, end, description=''):
    return dict(type=kind, description=description, location=dict(start=dict(value=start), end=dict(value=end)))


def entry(accession, gene, length, features=(), comments=()):
    return dict(primaryAccession=accession, uniProtkbId=accession + '_MOUSE', entryType='UniProtKB reviewed (Swiss-Prot)',
                genes=[dict(geneName=dict(value=gene))], organism=dict(taxonId=10090),
                proteinDescription=dict(recommendedName=dict(fullName=dict(value=gene + ' protein'))),
                entryAudit=dict(entryVersion=1, sequenceVersion=1), sequence=dict(value='A' * length),
                features=list(features), comments=list(comments))


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
}


def fake_fetch(url):
    gene = unquote(url).split('gene_exact:')[1].split(')')[0]
    return dict(results=ENTRIES.get(gene, [])), {'x-uniprot-release': '2026_05'}


class Inventory(unittest.TestCase):
    def test_flags_candidates_and_known_example_first(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root/'in').mkdir()
            pd.DataFrame(dict(ligand=['Lig1', 'Nrxn1', 'Lig1'], receptor=['Rec1', 'Rec2', 'CX'],
                              type=['diffusion', 'contact', 'diffusion'])).to_csv(root/'in'/'lr.csv', index=False)
            pd.DataFrame(dict(entity=['Lig1', 'Nrxn1', 'Rec1', 'Rec2', 'CX'],
                              subunits=['Lig1', 'Nrxn1', 'Rec1', 'Rec2', 'Rec1;Sub2'])).to_csv(root/'in'/'entities.csv',
                                                                                                 index=False)
            (root/'in'/'preparation.json').write_text(json.dumps(dict(sample='synthetic')))
            inv.run(argparse.Namespace(input=str(root/'in'), output=str(root/'out'), workers=2), fetch=fake_fetch)
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
            self.assertEqual(set(candidates.gene), {'Nrxn1', 'Rec1', 'Lig1'})
            meta = json.loads((root/'out'/'inventory.json').read_text())
            self.assertEqual(meta['uniprot_releases'], ['2026_05'])
            self.assertEqual(meta['genes_with_several_entries'], ['Nrxn1'])


if __name__ == '__main__':
    unittest.main()
