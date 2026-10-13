"""Scientific Y regressions: small controls, not a production sweep."""
import sys,math,json,tempfile,unittest
from pathlib import Path
from dataclasses import replace
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import second_theorem_Y as y

class YScience(unittest.TestCase):
    def test_linton_reference(self):
        c=y.lattice.lattice_sums(1.,2.,np.sqrt(2.),200,20)
        self.assertLess(abs(y.lattice.greens_periodic(0.,.01,c,2.,1.)-complex(-.4595298795,-.3509130869)),8e-11)

    def test_batched_matches_scalar(self):
        for k,e,a in [(2.8+.01j,.1,-.02),(3.14155,.02,-.0046),(2.9,.2,-.046)]:
            c=replace(y.CONFIG,bem_M=8)
            A=y.assemble_full_matrix(k,e,a,8,c);B=y.assemble_full_matrix_scalar(k,e,a,8,c)
            np.testing.assert_allclose(A,B,rtol=3e-9,atol=1e-9)

    def test_green_modal_harmonics_and_periodicity(self):
        for k in (2.8,3.14155,2.8+.01j):
            vals=[]
            for H in (20,40):
                g,r=y.make_green_functions(k,replace(y.CONFIG,harmonic_order=H))
                X=.35;yy=.8;eta=.75
                n=np.arange(1,2000);p=n*np.pi/2
                q=np.where(p<k.real,-1j*np.sqrt(k*k-p*p+0j),np.sqrt(p*p-k*k+0j))
                reference=-np.sum(np.exp(-q*X)*np.sin(p*(yy+1))*np.sin(p*(eta+1))/q)/2
                self.assertLess(abs(g(X,yy,0,eta)-reference)/max(abs(reference),1),1e-10)
                for kernel in (g,r):
                    self.assertLess(abs(kernel(X,yy,0,eta)-kernel(X,yy+4,0,eta-4)),2e-11)
                vals.append(g(X,yy,0,eta))
            self.assertLess(abs(vals[0]-vals[1]),1e-11*max(abs(vals[1]),1))

    def test_modal_operator_independent(self):
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'runners'))
        from check_X_modal_operator import modal_matrix
        c=replace(y.CONFIG,bem_M=8);a=-.02;e=.1;k=2.8
        ref=modal_matrix(k,e,c,12000,geometry=lambda t,e,c:y.obstacle_geometry(t,e,a,c))
        A=y.assemble_full_matrix(k,e,a,8,c)
        self.assertLess(np.linalg.norm(A-ref)/np.linalg.norm(ref),1e-8)

    def test_offgrid_batched_matches_scalar(self):
        e=.1;a=-.02;k=2.8;c=y.CONFIG
        source=y.boundary_nodes(8);target=source+.07
        xx,yy,*_=y.obstacle_geometry(target,e,a,c)
        xi,eta,xp,yp,*_=y.obstacle_geometry(source,e,a,c)
        g,r=y.make_green_functions(k,c)
        actual=y._offgrid_kernel_array(xx,yy,xi,eta,xp,yp,k,c)
        reference=np.array([[y.weighted_normal_kernel(float(t),float(u),e,a,c,g,r) for u in source] for t in target])
        np.testing.assert_allclose(actual,reference,rtol=1e-8,atol=1e-9)

    def test_false_scale_and_nan_certificate(self):
        A=np.eye(24);A[0,0]=1e10
        with patch.object(y,'assemble_full_matrix',return_value=A):
            self.assertFalse(y.certify_full_mode(2.8,.1,-.02,24,True,1e4,y.CONFIG)['supported'])
        with patch.object(y,'assemble_full_matrix',return_value=np.full((24,24),np.nan)):
            with self.assertRaises(ValueError):y.certify_full_mode(2.8,.1,-.02,24,True,1e4,y.CONFIG)

    def test_nan_and_failure_do_not_approve_scan(self):
        mode=type('Mode',(),{'a_numerical':-.02,'kb_numerical':3.,'normalized_residual':1e-12})()
        for effect in (None,ValueError('injected')):
            with patch.object(y,'build_uniqueness_grid',return_value=np.linspace(1.8,2.8,9)),patch.object(y,'full_singular_metrics',return_value=(math.nan,1.,math.nan,math.nan),side_effect=effect):
                result=y.whole_band_uniqueness_screen(.1,1.,mode,y.CONFIG)
            self.assertFalse(result.diagnostics['complete'])
            self.assertFalse(result.diagnostics['sampled_screen_passed'])
            self.assertFalse(result.diagnostics['uniqueness_supported'])

    def test_candidate_limit_not_silently_discarded(self):
        mode=type('Mode',(),{'a_numerical':-.02,'kb_numerical':3.,'normalized_residual':1e-12})()
        with patch.object(y,'build_uniqueness_grid',return_value=np.linspace(1.8,2.8,9)),patch.object(y,'full_singular_metrics',return_value=(1.,1.,1.,1.)),patch.object(y,'sampled_local_minimum_brackets',return_value=[(1.9,2.1),(2.2,2.4)]),patch.object(y,'refine_k_bracket',side_effect=RuntimeError('failure')):
            result=y.whole_band_uniqueness_screen(.1,1.,mode,replace(y.CONFIG,uniqueness_max_candidates=1))
        self.assertTrue(any(r.status=='pending/candidate-limit' for r in result.candidates))
        self.assertFalse(result.diagnostics['sampled_screen_passed'])

    def test_optimizer_failure_not_accepted(self):
        fake=type('R',(),{'success':False,'x':np.zeros(2),'message':'limit','nfev':80})()
        with patch.object(y,'least_squares',return_value=fake),patch.object(y,'near_zero_eigenpair',return_value=(1e-12+1e-12j,np.ones(24)/np.sqrt(24),1.)),patch.object(y,'_root_residual',return_value=np.zeros(2)),patch.object(y,'_root_jacobian',return_value=np.eye(2)):
            result,*_=y._solve_complex_eigenvalue_root(.1,-.2,.2,2.,[('leading',np.zeros(2))],y.CONFIG)
        self.assertFalse(result.success)
        self.assertFalse(result.attempts[0]['accepted'])

    def test_refinement_failure_rejected(self):
        fake=type('R',(),{'success':False,'message':'limit','x':.1})()
        with patch.object(y,'full_singular_metrics',return_value=(1.,1.,1.,1.)),patch.object(y,'minimize_scalar',return_value=fake):
            with self.assertRaises(RuntimeError):y.refine_k_bracket(2.,2.1,.1,-.02,y.CONFIG)

    def test_no_silent_sigma_clipping(self):
        with self.assertRaises(ValueError):y._joint_candidate_from_variables(.1,-.2,100.,0.,0.,y.CONFIG)
        with patch.object(y,'near_zero_eigenpair',side_effect=ValueError('injected')):
            with self.assertRaises(ValueError):y._root_residual(np.zeros(2),.1,-.2,.2,y.CONFIG)

    def test_mfs_circle_and_small_beta(self):
        circle,*_=y._solve_reference_mfs(y.CONFIG,80,beta_override=0.)
        self.assertLess(abs(circle.mu-1),1e-8);self.assertLess(abs(circle.a1_formula_2_7),1e-8)
        for beta in (.001,.02,.5):
            row,*_=y._solve_reference_mfs(y.CONFIG,120,beta_override=beta)
            self.assertLess(row.offgrid_relative_residual,1e-7)
            self.assertLess(row.identity_relative_error,1e-7)
            if beta<.05:self.assertLess(abs(row.a1_formula_2_7/beta+.5),1e-4)

    def test_printed_example_integral(self):
        t=2*np.pi*(np.arange(4096)+.5)/4096
        # Exact beta coefficient of the printed numerator: N1.
        def integral(b):
            X=np.sin(t)-b/2*np.sin(2*t);Y=-np.cos(t)+b/2*np.cos(2*t)
            xp=np.cos(t)-b*np.cos(2*t);yp=np.sin(t)-b*np.sin(2*t)
            f=-np.cos(t)+b*(.5*np.cos(2*t)-.25)
            return np.mean((Y*xp-3*X*yp)*f)*2/3
        h=1e-5
        self.assertAlmostEqual((integral(h)-integral(-h))/(2*h),-.5,places=9)

    def test_signed_channel_recovers_formula_2_7(self):
        # Leading boundary density L0(a1+Y)=a1/2+(Y-Psi)/2.
        # Expanding its EXACT outgoing coefficient recovers (5.47)/(2.7).
        row,c,sx,sy=y._solve_reference_mfs(y.CONFIG,120)
        t=y.boundary_nodes(1024);X,Y,*_=y.reference_obstacle_geometry(t,y.CONFIG)
        psi=-np.log(np.hypot(X[:,None]-sx,Y[:,None]-sy))@c/(2*np.pi)
        density=.5*(Y-psi+row.a1_formula_2_7)
        errors=[]
        for e in (.02,.01,.005):
            _,coef=y.modal_open_channel_coefficients(np.pi,e,e*row.a1_formula_2_7,density.astype(complex),y.CONFIG)
            q=np.sqrt(np.pi**2-y.lambda_1(y.CONFIG))
            errors.append(abs(1j*2*q*coef/e**2))
        self.assertLess(errors[-1],2e-5)
        self.assertLess(errors[1]/errors[0],.26)
        self.assertLess(errors[2]/errors[1],.26)

    def test_direct_open_coefficient_and_modal_field(self):
        k=2.8;e=.1;a=-.02;t=y.boundary_nodes(24);v=np.cos(t).astype(complex)
        left,right=y.modal_open_channel_coefficients(k,e,a,v,y.CONFIG)
        q=np.sqrt(k*k-y.lambda_1(y.CONFIG))
        nodes,weights=np.polynomial.legendre.leggauss(96);phi=np.cos(np.pi*nodes/2)
        g=y.make_extended_field_green(k,y.CONFIG)
        for xx,coef in [(-2.,left),(2.,right)]:
            field=y.reconstruct_field_at_points([(xx,float(yy)) for yy in nodes],v,t,e,a,y.CONFIG,g)
            projected=np.dot(weights*phi,field)/np.dot(weights,phi*phi)*np.exp(-1j*q*abs(xx))
            self.assertLess(abs(projected-coef),1e-12)
        # Analytical modal source derivative agrees with the independent FD route.
        point=(2.,.3)
        explicit=2*np.pi/24*sum(y.weighted_kernel_at_field_point(*point,float(tt),e,a,y.CONFIG,g)*vv for tt,vv in zip(t,v))
        actual=y.reconstruct_field_at_points([point],v,t,e,a,y.CONFIG,g)[0]
        self.assertLess(abs(explicit-actual),2e-10)

    def test_missing_physics_never_supported(self):
        self.assertFalse(y.physical_is_supported(None))

    def test_primary_dedup_uses_local_resolution_and_vector(self):
        vector=np.ones(24)/np.sqrt(24)
        mode=type('Mode',(),{'a_numerical':-.023,'kb_numerical':3.1254703839625,'normalized_residual':3e-12})()
        with patch.object(y,'full_singular_pair',return_value=(0.,1.,0.,0.,vector)):
            self.assertTrue(y.candidate_matches_primary(mode.kb_numerical+6e-11,.318,vector,.1,mode,y.CONFIG))
            self.assertFalse(y.candidate_matches_primary(mode.kb_numerical+1e-5,.318,vector,.1,mode,y.CONFIG))
            orthogonal=np.arange(24)-11.5
            self.assertFalse(y.candidate_matches_primary(mode.kb_numerical+6e-11,.318,orthogonal,.1,mode,y.CONFIG))

if __name__=='__main__':unittest.main()
