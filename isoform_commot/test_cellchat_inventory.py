import unittest

import pandas as pd

from cellchat_inventory import complex_table, conflicts, sharing


def interactions(rows):
    it = pd.DataFrame(rows, columns=['ligand', 'receptor', 'annotation', 'pathway_name', 'evidence'])
    complexes = {'Cx'}
    it['ligand_complex'] = it.ligand.isin(complexes)
    it['receptor_complex'] = it.receptor.isin(complexes)
    return it


class InventoryTests(unittest.TestCase):
    def test_sharing_counts(self):
        p = pd.DataFrame(dict(ligand=['A', 'A', 'B', 'C'], receptor=['X', 'Y', 'Y', 'Z']))
        s = sharing(p)
        self.assertEqual((s['pairs'], s['ligands'], s['receptors']), (4, 3, 3))
        self.assertEqual(s['pairs_sharing_ligand'], 2)      # A-X, A-Y
        self.assertEqual(s['pairs_sharing_receptor'], 2)    # A-Y, B-Y
        self.assertEqual(s['pairs_sharing_either'], 3)      # C-Z shares nothing
        self.assertEqual(sharing(p.iloc[:0])['pairs'], 0)

    def test_complexes_counted_separately(self):
        it = interactions([['A', 'X', 'Secreted Signaling', 'P', ''], ['A', 'Cx', 'Secreted Signaling', 'P', ''],
                           ['Cx', 'X', 'Cell-Cell Contact', 'P', '']])
        t = complex_table(it).set_index('cls')
        self.assertEqual(t.loc['Secreted Signaling', 'single_gene'], 1)
        self.assertEqual(t.loc['Secreted Signaling', 'complex_receptor_only'], 1)
        self.assertEqual(t.loc['Cell-Cell Contact', 'complex_ligand_only'], 1)

    def test_conflicts_orientation_and_homophilic(self):
        it = interactions([['Nrxn1', 'Nlgn1', 'Cell-Cell Contact', 'NRXN', ''],
                           ['Negr1', 'Negr1', 'Cell-Cell Contact', 'NEGR', ''],
                           ['Vegfa', 'Kdr', 'Secreted Signaling', 'VEGF', ''],
                           ['Col1a1', 'Cd44', 'ECM-Receptor', 'COLLAGEN', '']])
        p = it[it.annotation != 'ECM-Receptor'].assign(step_b=True, step_c=False)
        ccc = pd.DataFrame(dict(ligand=['NLGN1', 'NEGR1', 'VEGFA', 'COL1A1'],
                                receptor=['NRXN1', 'NEGR1', 'KDR', 'CD44'],
                                type=['diffusion', 'diffusion', 'diffusion', 'diffusion'], origin='x'))
        m = conflicts(it, p, ccc).set_index(['ligand', 'receptor'])
        self.assertEqual(m.loc[('Nrxn1', 'Nlgn1'), 'orientation'], 'swapped')
        self.assertEqual(m.loc[('Nrxn1', 'Nlgn1'), 'relation'], 'conflict')
        self.assertEqual(len(m.loc[[('Negr1', 'Negr1')]]), 1)  # not duplicated as "swapped"
        self.assertEqual(m.loc[('Vegfa', 'Kdr'), 'relation'], 'agree')
        self.assertEqual(m.loc[('Col1a1', 'Cd44'), 'relation'], 'cellchat_class_without_cccvelo_counterpart')


if __name__ == '__main__':
    unittest.main()
