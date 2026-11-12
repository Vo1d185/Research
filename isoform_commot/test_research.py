import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize, LinearConstraint
from research import solve, kernel, cost, demo, run


class ResearchTests(unittest.TestCase):
    def test_kernel_units_and_limit(self):
        d = np.array([1., 2., 4.])
        for family in ['inverse','r0_raw','r0_shape']:
            np.testing.assert_allclose(kernel(d, family, 1, 2),
                                       kernel(d*1000, family, 1000, 2000))
        np.testing.assert_allclose(kernel(d,'r0_shape',1,1e-9), kernel(d,'inverse',1), atol=1e-9)

    def test_contact_zero_cost_and_exclusions(self):
        D = np.array([[0.,1.,2.],[1.,0.,1.],[2.,1.,0.]])
        pairs = pd.DataFrame([['L','R','contact']], columns=['ligand','receptor','type'])
        C = cost(D,pairs,['L'],['R'],[1.], 'inverse',None,1.,1.1,5.)
        self.assertEqual(C[0,1],0)
        self.assertTrue(np.isinf(C[0,0]) and np.isinf(C[0,2]))
        for invalid in [np.nan, -1, 1.1]:
            with self.assertRaises(ValueError):
                cost(D,pairs,['L'],['R'],[invalid], 'inverse',None,1.,1.1,5.)

    def test_empty_support(self):
        np.testing.assert_array_equal(solve(np.full((2,2),np.inf), [1,2],[2,1]), np.zeros((2,2)))

    def test_ppi_changes_competition(self):
        C = np.array([[0.,.5],[.5,0.]])
        before = solve(C,[.5,.5],[.5,.5])
        C[0,0] -= np.log(.1)
        after = solve(C,[.5,.5],[.5,.5])
        self.assertLess(after[0,0],before[0,0]-1e-6)

    def test_primal_oracle(self):
        a,b = np.array([.45,.55]), np.array([.6,.4])
        C = np.array([[.2,.6],[.7,.1]])
        eps,rho = .5,.7
        T = solve(C,a,b,eps,rho)
        def entropy(z):
            z=np.maximum(z,1e-15)
            return (z*(np.log(z)-1)).sum()
        def objective(z):
            t=z.reshape(2,2); mu=a-t.sum(1); nu=b-t.sum(0)
            return (C*t).sum()+rho*(mu.sum()+nu.sum())+eps*(entropy(t)+entropy(mu)+entropy(nu))
        A=np.array([[1,1,0,0],[0,0,1,1],[1,0,1,0],[0,1,0,1]])
        opt=minimize(objective, np.full(4,.05), method='SLSQP', bounds=[(1e-12,None)]*4,
                     constraints=LinearConstraint(A,-np.inf,np.r_[a,b]-1e-12),
                     options={'ftol':1e-12,'maxiter':1000})
        self.assertTrue(opt.success,opt.message)
        np.testing.assert_allclose(T,opt.x.reshape(2,2),atol=1e-4)

    def test_pipeline_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            src,out=Path(d)/'input',Path(d)/'output'
            demo(src); run(src,out)
            summary=pd.read_csv(out/'summary.csv')
            self.assertEqual(set(summary.family), {'inverse','r0_raw','r0_shape'})
            self.assertIn('ppi',set(summary.score_mode))
            self.assertTrue((summary.total_T >= 0).all())
            with self.assertRaises(FileExistsError):
                demo(src)
            with self.assertRaises(FileExistsError):
                run(src,out)
            ppi=pd.read_csv(src/'ppi.csv').iloc[:-1]
            ppi.to_csv(src/'ppi.csv',index=False)
            with self.assertRaisesRegex(ValueError,'coverage'):
                run(src,Path(d)/'missing')


if __name__ == '__main__':
    unittest.main()
