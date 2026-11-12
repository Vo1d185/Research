import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from analyze import paired_metrics, type_and_pair_metrics, analyze
from research import demo, run


class AnalysisTests(unittest.TestCase):
    def test_scaling_is_not_reallocation(self):
        base = np.array([[1.,2.],[3.,4.]])
        result = paired_metrics(base,base/2)
        self.assertEqual(result['total_change_percent'],-50)
        self.assertEqual(result['paired_relative_L1'],.5)
        self.assertEqual(result['allocation_TV'],0)

    def test_disjoint_allocation_and_zero_mass(self):
        result = paired_metrics(np.array([1.,0.]), np.array([0.,1.]))
        self.assertEqual(result['total_change_percent'],0)
        self.assertEqual(result['paired_relative_L1'],2)
        self.assertEqual(result['allocation_TV'],1)
        result = paired_metrics(np.zeros(2), np.ones(2))
        self.assertIsNone(result['allocation_TV'])
        self.assertIsNone(result['paired_relative_L1'])

    def test_type_partition_and_distances(self):
        T = np.zeros((4,2))
        T[0,1] = 3
        T[3,0] = 2
        pairs = pd.DataFrame([['L1','R','diffusion'],['L2','R','contact']],
                             columns=['ligand','receptor','type'])
        rows, _ = type_and_pair_metrics(T,np.zeros_like(T),['a','b'],
            ['L1','L2'],['R'],pairs,np.array([[0.,2.],[2.,0.]]))
        self.assertEqual(rows[0]['total_T'],3)
        self.assertEqual(rows[1]['total_T'],2)
        self.assertEqual(rows[0]['mean_distance'],2)
        self.assertEqual(rows[1]['mean_distance'],2)

    def test_complete_analysis_and_input_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src,res,out = root/'input',root/'results',root/'analysis'
            demo(src); run(src,res); analyze(res,src,out)
            pairs = pd.read_csv(out/'ppi_matched.csv')
            self.assertEqual(len(pairs),9)
            summary = pd.read_csv(res/'summary.csv').set_index('case')
            for row in pairs.itertuples():
                self.assertEqual(summary.loc[row.baseline_case,'score_mode'],'no_ppi')
                self.assertEqual(summary.loc[row.baseline_case,'family'],row.family)
                np.testing.assert_allclose(summary.loc[row.baseline_case,'r0'],row.r0,equal_nan=True)
            typed = pd.read_csv(out/'by_type.csv').groupby('case').total_T.sum()
            np.testing.assert_allclose(typed.loc[summary.index],summary.total_T)
            self.assertEqual(json.loads((out/'analysis_manifest.json').read_text())['status'],'complete')
            for filename in ['distance_by_type.png','ppi_matched.png','lr_pair_changes.png']:
                self.assertGreater((out/filename).stat().st_size,1000)
            with self.assertRaises(FileExistsError):
                analyze(res,src,out)
            with (src/'spatial.csv').open('a') as f:
                f.write('\n')
            with self.assertRaisesRegex(ValueError,'Input changed'):
                analyze(res,src,root/'other_analysis')

    def test_no_ppi_input_does_not_invent_comparison(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            demo(root/'input')
            (root/'input'/'ppi.csv').unlink()
            run(root/'input',root/'results')
            analyze(root/'results',root/'input',root/'analysis')
            self.assertTrue(pd.read_csv(root/'analysis'/'ppi_matched.csv').empty)


if __name__ == '__main__':
    unittest.main()
