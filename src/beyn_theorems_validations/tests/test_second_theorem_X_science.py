"""Scientific regressions for X; no production sweep."""
import sys,math,tempfile,unittest
from pathlib import Path
from dataclasses import replace
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import second_theorem_X as x

class XScience(unittest.TestCase):
    def test_batched_matches_scalar_kernel(self):
        c=x.CONFIG;M=8;e=.1;k=2.8+.01j;t=x.boundary_nodes(M)
        g,r=x.make_green_functions(k,c)
        ref=np.eye(M,dtype=complex)
        for i,a in enumerate(t):
            for j,b in enumerate(t):ref[i,j]-=4*np.pi/M*x.weighted_normal_kernel(a,b,e,c,g,r)
        np.testing.assert_allclose(x.assemble_full_matrix(k,e,M,c),ref,rtol=1e-10,atol=1e-10)

    def test_green_modal_and_periodicity(self):
        g,_=x.make_green_functions(2.8,x.CONFIG)
        n=np.arange(1,1000);p=n*np.pi/2;k=2.8
        q=np.where(p<k,-1j*np.sqrt(k*k-p*p+0j),np.sqrt(p*p-k*k+0j))
        X=.5;y=.8;eta=.75
        ref=-np.sum(np.exp(-q*X)*np.sin(p*(y+1))*np.sin(p*(eta+1))/q)/2
        self.assertLess(abs(g(X,y,0,eta)-ref),1e-9)
        self.assertLess(abs(g(X,y,0,eta)-g(X,y+4,0,eta)),1e-13)

    def test_independent_accelerated_modal_matrix(self):
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'runners'))
        from check_X_modal_operator import modal_matrix
        c=replace(x.CONFIG,bem_M=8)
        A=x.assemble_odd_matrix(2.8,.1,8,c)
        B=modal_matrix(2.8,.1,c,12000)
        self.assertLess(np.linalg.norm(A-B)/np.linalg.norm(B),1e-8)

    def test_projection_and_full_mode_residual(self):
        Q=x.odd_lift_matrix(24);R=x.reflection_matrix(24)
        np.testing.assert_allclose(Q.conj().T@Q,np.eye(12),atol=1e-15)
        np.testing.assert_allclose(R@Q,-Q,atol=1e-15)
        B=np.eye(12);B[0]=B[1];B[0,0]=1e-12
        A=Q@B@Q.conj().T+(np.eye(24)-Q@Q.conj().T)
        with patch.object(x,'assemble_full_matrix',return_value=A):
            cert=x.certify_odd_mode(2.,.1,24,True,1e4,x.CONFIG)
        self.assertTrue(cert['supported'])

    def test_scale_disguised_valley_rejected(self):
        Q=x.odd_lift_matrix(24);B=np.eye(12);B[0,0]=1e10
        A=Q@B@Q.conj().T+(np.eye(24)-Q@Q.conj().T)
        with patch.object(x,'assemble_full_matrix',return_value=A):
            self.assertFalse(x.certify_odd_mode(2.,.1,24,True,1e4,x.CONFIG)['supported'])

    def test_nan_matrix_rejected(self):
        with patch.object(x,'assemble_full_matrix',return_value=np.full((24,24),np.nan)):
            with self.assertRaises(ValueError):x.certify_odd_mode(2.,.1,24,True,1e4,x.CONFIG)

    def test_failed_optimizer_rejected(self):
        fake=type('R',(),{'success':False,'x':.02,'message':'injected'})()
        with patch.object(x,'absolute_sv_from_sigma',return_value=1.),patch.object(x,'minimize_scalar',return_value=fake):
            with self.assertRaises(RuntimeError):x.refine_candidate_for_M(.1,24,.01,.03,x.CONFIG)

    def test_mu_circle_and_identity(self):
        row=x.compute_mu_mfs(x.CONFIG,80,beta_override=0.)
        self.assertLess(abs(row.mu-1),1e-8)
        self.assertLess(row.offgrid_relative_residual,1e-7)
        self.assertLess(row.identity_relative_error,1e-8)

    def test_modal_far_field(self):
        g=x.make_extended_field_green(2.8,x.CONFIG)
        self.assertTrue(np.isfinite(g(8.,.8,-.2,.1)))
        self.assertLess(abs(g(8.,1.,-.2,.1)),1e-14)

    def test_odd_green_matches_half_channel_modal(self):
        g,_=x.make_green_functions(2.8,x.CONFIG)
        xx=.35;y=.22;eta=.17;k=2.8
        n=np.arange(1,1000);p=n*np.pi;q=np.sqrt(p*p-k*k+0j)
        half=-np.sum(np.exp(-q*xx)*np.sin(p*y)*np.sin(p*eta)/q)
        self.assertLess(abs(g(xx,y,0,eta)-g(xx,y,0,-eta)-half),1e-9)

    def test_physical_failure_and_rank_disagreement_block_uniqueness(self):
        from dataclasses import fields
        diag=x.BeynDiagnostics(.1,24,12,384,1,2,1,1e9,1.,1e-14,1,1,0)
        diag.moments_converged=True;diag.candidate_stable=True;diag.rank_stable=True
        conv=[x.BeynConvergenceRow(.1,N,1,1.,1e-10,1e-12,1e-5,1,2.,0.,0.,True) for N in (192,384)]
        mode=x.ModeResult(.1,1,'mock',2.,2.,1e-12,1.,1e-12,1e4,math.nan,True)
        physical=x.PhysicalDiagnostics(**{f.name:False if f.type=='bool' else .1 for f in fields(x.PhysicalDiagnostics)})
        with tempfile.TemporaryDirectory() as tmp,patch.object(x,'run_beyn_convergence_study',return_value=(np.array([2.]),np.array([1.]),diag,[],[x.DiscoverySeed(2.,0.,'mock')],[diag],[],conv,math.nan)),patch.object(x,'run_candidate_refinement',return_value=([],mode)),patch.object(x,'beyn_discover',return_value=(np.array([2.]),np.array([1.]),diag,[],np.eye(12))),patch.object(x,'parity_commutator_residual',return_value=0.),patch.object(x,'physical_mode_diagnostics',return_value=physical):
            result,*_=x.validate_epsilon(.1,1.,0.,replace(x.CONFIG,output_directory=tmp))
        self.assertFalse(result.mode_exists_supported)
        self.assertFalse(result.odd_uniqueness_supported)
        self.assertFalse(result.unique_embedded_mode_supported)
        self.assertFalse(result.full_problem_uniqueness_supported)

    def test_disappearing_strict_candidate_is_preserved(self):
        def discover(e,N,mu,c,**kwargs):
            d=x.BeynDiagnostics(e,24,12,N,0 if N==64 else 1,0 if N==64 else 1,0 if N==64 else 1,1e9,1.,1e-15,1,1,0)
            d._s0=np.full((12,8),float(N),complex);d._s1=2*d._s0
            if N==64:return np.array([]),np.array([1.]),d,[],d._s0
            row=x.BeynEigenvalueRow(e,N,0,2.,0.,True,True,True,0.,False,True)
            return np.array([2.]),np.array([1.]),d,[row],d._s0
        c=replace(x.CONFIG,beyn_quadrature_levels=(32,64),beyn_adaptive_quadrature_levels=())
        with patch.object(x,'beyn_discover',side_effect=discover):
            result=x.run_beyn_convergence_study(.1,1.,c)
        self.assertEqual(len(result[4]),1)
        self.assertEqual(result[4][0].real,2.)

    def test_reduced_complex_conjugation(self):
        plus=x.assemble_odd_matrix(2.8+.01j,.1,24,x.CONFIG)
        minus=x.assemble_odd_matrix(2.8-.01j,.1,24,x.CONFIG)
        self.assertLess(np.linalg.norm(plus.conj()-minus)/np.linalg.norm(plus),1e-8)

    def test_graded_quadrature_counts_zero_one_two(self):
        angles,weights=x.contour_quadrature(512)
        self.assertAlmostEqual(float(np.sum(weights)),2*np.pi,places=13)
        self.assertTrue(np.all(weights>0))
        for roots in ((),(2.,),(2.,2.7)):
            def matrix(z,*args):
                A=np.eye(12,dtype=complex)
                for i,r in enumerate(roots):A[i,i]=z-r
                return A
            with patch.object(x,'assemble_odd_matrix',side_effect=matrix):
                eigs,_,diag,_,_=x.beyn_discover(.1,512,1.,x.CONFIG)
            self.assertEqual(diag.estimated_rank,len(roots))
            if roots:np.testing.assert_allclose(sorted(eigs.real),roots,rtol=1e-10,atol=1e-10)

    def test_moments_missing_do_not_converge(self):
        # Candidate+rank agreement alone must trigger adaptive levels.
        seen=[]
        def discover(e,N,mu,c,**kwargs):
            seen.append(N)
            d=x.BeynDiagnostics(e,24,12,N,1,1,1,1e9,1.,1e-15,1,1,0)
            d._s0=np.full((12,8),float(N),complex);d._s1=2*d._s0
            row=x.BeynEigenvalueRow(e,N,0,2.,0.,True,True,True,0.,False,True)
            return np.array([2.]),np.array([1.]),d,[row],d._s0
        c=replace(x.CONFIG,beyn_quadrature_levels=(32,64),beyn_adaptive_quadrature_levels=(128,))
        with patch.object(x,'beyn_discover',side_effect=discover):
            result=x.run_beyn_convergence_study(.1,1.,c)
        self.assertEqual(seen,[32,64,128]);self.assertFalse(result[2].moments_converged)

if __name__=='__main__':unittest.main()
