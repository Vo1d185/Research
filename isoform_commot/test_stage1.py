"""Synthetic fixtures test data handling; these are not biological results."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist

import stage1
from research import contact_neighbors, cost


class StageOneTests(unittest.TestCase):
    def fixture(self, root, expression_state='counts'):
        xy = np.array([[0, 0], [1, 0], [0, 1], [1.1, 1.2], [2, .5], [2, 2]], float)
        genes = ['Cxcl12', 'Cxcr4', 'Cxcr3', 'Dll1', 'Notch1', 'Notch2', 'Other']
        x = np.arange(1, 43).reshape(6, 7).astype(float)
        data = ad.AnnData(sparse.csr_matrix(x), obs=pd.DataFrame(index=[f's{i}' for i in range(6)]),
                         var=pd.DataFrame(index=genes))
        data.obsm['spatial'] = xy
        data.obs['scc_anno'] = ['A', 'B'] * 3
        data.write_h5ad(root/'fixture.h5ad')
        definitions = {
            'diff_LigRecDB.csv': ('diffusion', [('CXCL12', 'CXCR4'), ('CXCL12', 'CXCR3')]),
            'cont_LigRecDB.csv': ('contact', [('DLL1', 'NOTCH1'), ('DLL1', 'NOTCH2'),
                                            ('DLL1', 'Some_complex')]),
        }
        root.joinpath('lr').mkdir()
        checksums = {}
        for filename, (kind, rows) in definitions.items():
            path = root/'lr'/filename
            pd.DataFrame(rows, columns=['source', 'target']).assign(origin='test_fixture').to_csv(path, index=False)
            checksums[filename] = (kind, stage1.digest(path))
        args = stage1.parser().parse_args(['prepare', '--h5ad', str(root/'fixture.h5ad'),
            '--lr-dir', str(root/'lr'), '--output', str(root/'input'), '--matrix', 'X',
            '--expression-state', expression_state, '--symbol-policy', 'author-uppercase',
            '--coordinate-unit', 'synthetic', '--max-spots', '6', '--pairs-per-type', '2'])
        return x, xy, args, checksums

    def test_full_library_normalization_exact_lr_types_and_complete_pilot(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            x, _, args, checksums = self.fixture(root)
            with patch.dict(stage1.LR_FILES, checksums, clear=True), contextlib.redirect_stdout(io.StringIO()):
                stage1.prepare(args)
            got = pd.read_csv(root/'input'/'expression.csv', index_col=0)
            expected = np.log1p(x/x.sum(axis=1)[:, None]*10000)
            np.testing.assert_allclose(got['CXCL12'], expected[:, 0])
            self.assertNotIn('OTHER', got)
            pairs = pd.read_csv(root/'input'/'lr.csv')
            self.assertEqual(pairs.groupby('type').size().to_dict(), {'contact': 2, 'diffusion': 2})
            audit = pd.read_csv(root/'input'/'lr_audit.csv')
            self.assertEqual(audit.iloc[-1].reason, 'unsupported_complex_or_nonstandard_entity')
            report = json.loads((root/'input'/'preparation.json').read_text())
            self.assertEqual(report['shared_ligands'], 2)
            run_args = stage1.parser().parse_args(['run', '--input', str(root/'input'), '--output', str(root/'experiment')])
            with contextlib.redirect_stdout(io.StringIO()):
                stage1.run_stage(run_args)
            summary = pd.read_csv(root/'experiment'/'all_conditions.csv')
            self.assertEqual(len(summary), 5)
            self.assertEqual(set(summary.score_mode), {'no_ppi'})
            self.assertEqual(set(summary.family), {'inverse', 'r0_raw'})
            stability = pd.read_csv(root/'experiment'/'ot_parameter_stability.csv')
            self.assertTrue((stability.allocation_TV_vs_eps02_rho1 == 0).all())
            self.assertIsNone(json.loads((root/'experiment'/'stage1_manifest.json').read_text())['best_r0'])
            with (root/'input'/'lr.csv').open('a') as stream:
                stream.write('FAKE,PAIR,diffusion\n')
            run_args.output = str(root/'changed')
            with self.assertRaisesRegex(ValueError, 'Prepared input changed'):
                stage1.run_stage(run_args)
            self.assertFalse((root/'changed').exists())

    def test_processed_expression_is_not_normalized_again(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            x, _, args, checksums = self.fixture(root, 'processed')
            with patch.dict(stage1.LR_FILES, checksums, clear=True), contextlib.redirect_stdout(io.StringIO()):
                stage1.prepare(args)
            got = pd.read_csv(root/'input'/'expression.csv', index_col=0)
            np.testing.assert_array_equal(got['NOTCH1'], x[:, 4])
            with self.assertRaisesRegex(ValueError, 'Output exists'):
                stage1.prepare(args)

    def test_bad_mapping_and_unavailable_types_are_not_guessed(self):
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            stage1.map_genes(['Gene', 'GENE'], 'author-uppercase')
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _, _, args, checksums = self.fixture(root)
            args.symbol_policy = 'exact'
            with patch.dict(stage1.LR_FILES, checksums, clear=True):
                with self.assertRaisesRegex(ValueError, 'Insufficient LR type coverage'):
                    stage1.prepare(args)
            self.assertEqual(json.loads((root/'input'/'preparation.json').read_text())['status'], 'blocked')

    def test_refuse_false_counts_and_modified_lr_source(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _, _, args, checksums = self.fixture(root)
            data = ad.read_h5ad(args.h5ad); data.X.data += .3; data.write_h5ad(args.h5ad)
            with patch.dict(stage1.LR_FILES, checksums, clear=True):
                with self.assertRaisesRegex(ValueError, 'Non-integer'):
                    stage1.prepare(args)
                with (root/'lr'/'cont_LigRecDB.csv').open('a') as f:
                    f.write('FAKE,PAIR,test\n')
                with self.assertRaisesRegex(ValueError, 'Wrong/missing'):
                    stage1.load_lr(root/'lr')

    def test_contact_eligibility_is_preserved_in_cost(self):
        xy = np.array([[0, 0], [1, 0], [.4, .8], [1, 1], [8, 8]])
        D = cdist(xy, xy)
        edges = contact_neighbors(xy, D, 1.5, 'delaunay')
        self.assertFalse(edges[-1].any())
        self.assertFalse(np.diag(edges).any())
        np.testing.assert_array_equal(edges, edges.T)
        pairs = pd.DataFrame([['L', 'R', 'contact']], columns=['ligand', 'receptor', 'type'])
        C = cost(D, pairs, ['L'], ['R'], [1], 'r0_raw', 1, D[D > 0].min(), 1.5, np.inf, edges)
        np.testing.assert_array_equal(np.isfinite(C), edges)
        self.assertTrue((C[edges] == 0).all())
        with self.assertRaisesRegex(ValueError, 'non-collinear'):
            contact_neighbors(np.array([[0, 0], [1, 0], [2, 0]]), np.zeros((3, 3)), 1, 'delaunay')

    def test_download_rejects_html_and_hash_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            source = root/'source'; target = root/'target'
            source.write_text('<html>expired link</html>')
            with self.assertRaisesRegex(ValueError, 'HDF5'):
                stage1.download(source.as_uri(), target, hdf5=True)
            self.assertFalse(target.exists())
            self.assertFalse(target.with_suffix('.part').exists())
            with self.assertRaisesRegex(ValueError, 'checksum'):
                stage1.download(source.as_uri(), target, expected='wrong')

    def test_download_validates_complete_h5ad_before_publishing(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            original, target = root/'original.h5ad', root/'downloaded.h5ad'
            ad.AnnData(np.ones((3, 4))).write_h5ad(original)
            stage1.download(original.as_uri(), target, hdf5=True)
            self.assertEqual(stage1.digest(original), stage1.digest(target))
            # A truncated file still has the correct signature: reject it before
            # renaming, so a later inspect never sees a purported completed file.
            broken = root/'truncated.h5ad'
            broken.write_bytes(original.read_bytes()[:64])
            partial_target = root/'should_not_exist.h5ad'
            with self.assertRaisesRegex(ValueError, 'incomplete or unreadable'):
                stage1.download(broken.as_uri(), partial_target, hdf5=True)
            self.assertFalse(partial_target.exists())
            self.assertFalse(partial_target.with_name(partial_target.name+'.part').exists())


if __name__ == '__main__':
    unittest.main()
