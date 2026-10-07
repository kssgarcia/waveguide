"""Small numerical and fail-closed regressions; never execute the production sweep."""
from __future__ import annotations

import contextlib
import io
import json
import math
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import first_theorem_beyn as f


def modal_green(X, y, eta, k, width=2.0, terms=2000):
    """Independent separation of variables, Linton (4.9). X is nonzero."""
    n = np.arange(1, terms + 1)
    p = n * np.pi / width
    q = np.sqrt(p * p - complex(k) ** 2 + 0j)
    return -np.sum(np.exp(-q * abs(X)) * np.sin(p * y) * np.sin(p * eta) / q) / width


def diagnostic(**changes):
    base = dict(
        epsilon=.1, M=24, quadrature_points=96, probe_dim=8,
        threshold_rank=1, gap_rank=1, estimated_rank=1, selected_gap_ratio=1e9,
        leading_s0_singular_value=1., trailing_kept_s0_singular_value=1.,
        first_discarded_s0_singular_value=1e-12, max_linear_solve_relative_residual=1e-14,
        raw_eigenvalues=1, near_real_candidates=1, near_contour_seeds=0,
        rank_stable=True, candidate_stable=True, moments_converged=True,
    )
    return f.BeynDiagnostics(**(base | changes))


class ContourQuadratureTests(unittest.TestCase):
    def test_positive_weights_cover_circle_and_smooth_integral(self):
        t, w = f.contour_quadrature(1024, f.CONFIG)
        self.assertTrue(np.all(w > 0))
        self.assertAlmostEqual(float(w.sum()), 2 * np.pi, places=13)
        self.assertLess(abs(np.sum(w * np.exp(3j * t))), 1e-12)
        self.assertTrue(np.all((t > 0) & (t < 2 * np.pi)))

    def test_cutoff_near_pole_count_and_moment_converge(self):
        config = replace(f.CONFIG, beyn_cutoff_margin=1e-8)
        _, right, _, _, _ = f.contour_geometry(config)
        pole = right - 1e-7
        values = []
        for budget in (1024, 2048):
            t, w = f.contour_quadrature(budget, config)
            z = np.array([f.ellipse_point(float(v), config)[0] for v in t])
            dz = np.array([f.ellipse_point(float(v), config)[1] for v in t])
            m0 = np.sum(w * dz / (z - pole)) / (2j * np.pi)
            m1 = np.sum(w * dz * z / (z - pole)) / (2j * np.pi)
            values.append((m0, m1))
        self.assertLess(abs(values[-1][0] - 1), 1e-7)
        self.assertLess(abs(values[-1][1] - pole), 1e-7)
        self.assertLess(abs(values[-1][0] - values[0][0]), 1e-5)

    def test_midpoint_compatibility_and_invalid_rule(self):
        t, w = f.contour_quadrature(96, replace(f.CONFIG, beyn_quadrature_rule="periodic-midpoint"))
        self.assertEqual(len(t), 96)
        self.assertAlmostEqual(float(w.sum()), 2 * np.pi)
        with self.assertRaises(ValueError):
            f.contour_quadrature(96, replace(f.CONFIG, beyn_quadrature_rule="invalid"))


class GreenTests(unittest.TestCase):
    def test_linton_table_2(self):
        S = f.lattice.lattice_sums(1., 2., np.sqrt(2.), 200, 20)
        actual = f.lattice.greens_periodic(0., .01, S, 2., 1.)
        self.assertLess(abs(actual - complex(-.4595298795, -.3509130869)), 8e-11)

    def test_difficult_boundary_points_against_modal_and_harmonics(self):
        # Previously ~20% Green error in the near-wall arc with unwrapped images.
        for k in [1.56, 1.56+.002j, 1.56-.002j]:
            for a, e in [(.65, .297895), (.6, .348421), (.7, .247368)]:
                values = []
                for H in [20, 40]:
                    config = replace(f.CONFIG, a=a, harmonic_order=H)
                    G, _ = f.make_green_functions(k, config)
                    t, s = 1.1, 1.8
                    x, y = e*np.cos(t), e*np.sin(t)
                    xi, eta = e*np.cos(s), e*np.sin(s)
                    actual = G(x, y, xi, eta)
                    reference = modal_green(x-xi, y+1+a, eta+1+a, k)
                    self.assertLess(abs(actual-reference)/abs(reference), 2e-10)
                    values.append(actual)
                self.assertLess(abs(values[0]-values[1]), 1e-11*max(abs(values[1]), 1.))

    def test_periodic_and_regularized_are_consistent(self):
        config = f.CONFIG
        G, R = f.make_green_functions(1.2+.01j, config)
        for kernel in [G, R]:
            # Shift field +period and source -period together: direct separation
            # changes by 2 periods, image separation remains identical.
            self.assertLess(abs(kernel(.1,.05,-.07,-.03)-kernel(.1,4.05,-.07,-4.03)), 2e-12)
        x,y,xi,eta = .1,.05,-.07,-.03
        r = np.hypot(x-xi,y-eta)
        # Algebraic H0-J0 difference, valid also at complex k.
        expected = -.25j*(f.lattice.H0((1.2+.01j)*r)-f.lattice.J0((1.2+.01j)*r))
        self.assertLess(abs(G(x,y,xi,eta)-R(x,y,xi,eta)-expected), 2e-12)

    def test_vectorized_assembly_matches_scalar_midpoint_rule(self):
        for k in [1.5704, .8+.01j]:
            config = replace(f.CONFIG, bem_M=8)
            actual = f.assemble_matrix(k,.1,8,config)
            G,R = f.make_green_functions(k,config)
            theta = f.boundary_nodes(8)
            reference = np.eye(8,dtype=complex)
            for i,p in enumerate(theta):
                for j,q in enumerate(theta):
                    reference[i,j] -= 4*np.pi/8*f.weighted_normal_kernel(float(p),float(q),.1,config,G,R)
            np.testing.assert_allclose(actual,reference,rtol=3e-8,atol=3e-9)

    def test_field_and_assembly_use_same_green(self):
        G,_ = f.make_green_functions(1.2,f.CONFIG)
        field = f.make_extended_field_green(1.2,f.CONFIG)
        self.assertLess(abs(G(.2,.1,-.1,.03)-field(.2,.1,-.1,.03)), 1e-13)


class CertificateTests(unittest.TestCase):
    def test_large_unrelated_scale_does_not_hide_nonsingular_valley(self):
        with patch.object(f,'assemble_matrix',return_value=np.diag([1e10,1.]).astype(complex)):
            c = f.certify_local_mode(1.,.1,2,True,1e4,f.CONFIG)
        self.assertLess(c.normalized_residual,1e-8)
        self.assertGreater(c.equilibrated_residual,.5)
        self.assertFalse(c.supported)

    def test_small_ratio_but_wrong_vector_large_residual_rejected(self):
        A = np.diag([1e-12,1.]).astype(complex)
        # Mimic corrupted/stale singular-vector data; the explicit A@v must catch it.
        with patch.object(f,'assemble_matrix',return_value=A), patch.object(
            f.np.linalg,'svd',side_effect=[(np.eye(2),np.array([1.,1e-12]),np.eye(2)),np.ones(2)]
        ):
            c = f.certify_local_mode(1.,.1,2,True,1e4,f.CONFIG)
        self.assertGreater(c.normalized_residual,.5)
        self.assertFalse(c.supported)

    def test_genuine_near_null_vector_passes(self):
        A = np.array([[1.,1.],[1.,1.+1e-10]],complex)
        with patch.object(f,'assemble_matrix',return_value=A):
            self.assertTrue(f.certify_local_mode(1.,.1,2,True,1000.,f.CONFIG).supported)
            self.assertFalse(f.certify_local_mode(1.,.1,2,False,1000.,f.CONFIG).supported)
            self.assertFalse(f.certify_local_mode(1.,.1,2,True,1.,f.CONFIG).supported)

    def test_nonfinite_matrix_rejected(self):
        with patch.object(f,'assemble_matrix',return_value=np.full((2,2),np.nan)):
            self.assertFalse(f.certify_local_mode(1.,.1,2,True,1e3,f.CONFIG).supported)


class NonexistenceTests(unittest.TestCase):
    def scan(self, objective, **kwargs):
        with patch.object(f,'absolute_singular_value_from_sigma',objective), contextlib.redirect_stdout(io.StringIO()):
            return f.independent_whole_band_svd_scan(.1,.35,kwargs.get('config',f.CONFIG),'regression')[0]

    def test_all_nan_is_ambiguous(self):
        row = self.scan(lambda *args: math.nan)
        self.assertFalse(row.nonexistence_supported)
        self.assertFalse(row.diagnostic_complete)
        self.assertEqual(row.status,'ambiguous')

    def test_one_nan_is_ambiguous(self):
        calls = 0
        def objective(s,*args):
            nonlocal calls
            calls += 1
            return math.nan if calls==5 else 1+s
        self.assertFalse(self.scan(objective).nonexistence_supported)

    def test_evaluation_exception_is_ambiguous(self):
        def objective(*args): raise ValueError('injected evaluation failure')
        self.assertFalse(self.scan(objective).nonexistence_supported)

    def test_failed_refinement_is_ambiguous(self):
        with patch.object(f,'run_candidate_refinement',side_effect=RuntimeError('injected')):
            row=self.scan(lambda s,*args: 1+(math.log(s/.02))**2)
        self.assertGreater(row.failed_refinements,0)
        self.assertFalse(row.nonexistence_supported)

    def test_discarded_candidates_are_ambiguous(self):
        config=replace(f.CONFIG,nonexistence_scan_max_candidates=0)
        row=self.scan(lambda s,*args: 1+(math.log(s/.02))**2,config=config)
        self.assertGreater(row.unevaluated_candidates,0)
        self.assertFalse(row.nonexistence_supported)

    def test_unresolved_valley_is_ambiguous(self):
        mode=f.ModeResult(.1,1,'test',1.,0.,1.,.02,1.,1e10,1e-10,1e4,math.nan,False)
        with patch.object(f,'run_candidate_refinement',return_value=([],mode)):
            row=self.scan(lambda s,*args: 1+(math.log(s/.02))**2)
        self.assertFalse(row.nonexistence_supported)
        self.assertFalse(row.diagnostic_complete)

    def test_completed_monotone_scan_supports_non_detection(self):
        row=self.scan(lambda s,*args: 1+s)
        self.assertTrue(row.nonexistence_supported)
        self.assertTrue(row.diagnostic_complete)


class BeynTests(unittest.TestCase):
    def test_rank_alone_never_certifies_uniqueness(self):
        for overrides in [dict(moments_converged=False),dict(candidate_stable=False),dict(rank_stable=False),dict(threshold_rank=2)]:
            self.assertFalse(f.uniqueness_is_supported(diagnostic(**overrides),1,True,f.CONFIG))
        self.assertFalse(f.uniqueness_is_supported(diagnostic(),1,False,f.CONFIG))
        self.assertTrue(f.uniqueness_is_supported(diagnostic(),1,True,f.CONFIG))

    def test_both_moments_must_converge(self):
        self.assertFalse(f.moments_converged(1e-5,.1,f.CONFIG))
        self.assertFalse(f.moments_converged(math.nan,1e-5,f.CONFIG))
        self.assertTrue(f.moments_converged(1e-5,1e-5,f.CONFIG))

    def test_diagonal_nonlinear_eigenproblem(self):
        cfg=replace(f.CONFIG,bem_M=4,beyn_probe_dim=3)
        for roots in [[],[.5],[.5,.8]]:
            def assemble(z,e,M,c):
                values=np.ones(M,dtype=complex)
                for j,r in enumerate(roots): values[j]=z-r
                return np.diag(values)
            with patch.object(f,'assemble_matrix',assemble),contextlib.redirect_stdout(io.StringIO()):
                raw,_,diag,_,_,_=f.beyn_discover(.1,192,cfg)
            self.assertEqual(diag.estimated_rank,len(roots))
            np.testing.assert_allclose(np.sort(raw.real),roots,atol=1e-12)


class CacheTests(unittest.TestCase):
    def test_configuration_and_geometry_invalidate_fingerprint(self):
        base=f.configuration_fingerprint(f.CONFIG,.1,.6)
        for changes in [dict(bem_M=32),dict(harmonic_order=40),dict(lattice_terms=100),
                        dict(finite_difference_step=1e-7),dict(b=2.),dict(beyn_moment_convergence_tolerance=.01),
                        dict(singular_vector_residual_tolerance=1e-5),dict(paper_search_quadrature=192)]:
            self.assertNotEqual(base,f.configuration_fingerprint(replace(f.CONFIG,**changes),.1,.6))
        self.assertNotEqual(base,f.configuration_fingerprint(f.CONFIG,.2,.6))
        self.assertNotEqual(base,f.configuration_fingerprint(f.CONFIG,.1,.65))
        with patch.object(f,'ASSEMBLY_VERSION','future-version'):
            self.assertNotEqual(base,f.configuration_fingerprint(f.CONFIG,.1,.6))

    def test_legacy_directory_cannot_be_written(self):
        with self.assertRaises(ValueError): f.require_corrected_output(f.LEGACY_OUTPUT)
        with self.assertRaises(ValueError): f.require_corrected_output(f.LEGACY_OUTPUT/'new')

    def test_unsigned_legacy_cache_is_not_reused(self):
        with tempfile.TemporaryDirectory() as folder:
            directory=Path(folder)
            (directory/'paper_a_sweep_M24.csv').write_text('a,epsilon,M\n.6,.1,24\n')
            with patch.object(f,'run_paper_a_sweep',return_value=[]) as compute,contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(f.load_or_run_paper_a_sweep(f.CONFIG,directory),[])
            compute.assert_called_once()


class CrossRouteTests(unittest.TestCase):
    def test_local_existence_retained_when_global_moments_do_not_converge(self):
        config=replace(f.CONFIG,run_physical_diagnostics=False)
        def discover(e,c,margin,nq):
            diag=diagnostic(quadrature_points=nq)
            diag._s0=np.full((2,2),1.+nq/1000.,complex)
            diag._s1=2*diag._s0
            return c,np.array([1.]),diag,[],[complex(1.)]
        mode=f.ModeResult(.1,1,'test',1.,0.,1.,.02,1e-10,2.,5e-11,1e4,math.nan,True,1e-10,1e-10,1e-10)
        with patch.object(f,'_paper_discover_at_margin',discover),patch.object(
            f,'run_candidate_refinement',return_value=([],mode)
        ):
            row=f._paper_search_case(.1,.6,config,(1e-5,),True,1.,.02,2.)
        self.assertTrue(row.mode_exists_supported)
        self.assertFalse(row.uniqueness_supported)
        self.assertEqual(row.status,'mode-found')
        self.assertTrue(np.isfinite(row.kb))

    def test_candidate_lost_between_levels_cannot_support_absence(self):
        config=replace(f.CONFIG,run_physical_diagnostics=False)
        def discover(e,c,margin,nq):
            candidates=[complex(1.)] if nq==config.paper_verify_quadrature else []
            diag=diagnostic(estimated_rank=0,threshold_rank=0,quadrature_points=nq)
            diag._s0=np.zeros((2,2),complex);diag._s1=np.zeros((2,2),complex)
            return c,np.array(candidates),diag,[],candidates
        scan=f.IndependentScanSummary(.1,.6,'test',1e-6,1.,160,24,0,0,0,0,math.nan,math.nan,False,True,'no-resolved-interior-mode',True)
        with patch.object(f,'_paper_discover_at_margin',discover),patch.object(
            f,'candidate_sigma_tasks',return_value=[]
        ),patch.object(f,'independent_whole_band_svd_scan',return_value=(scan,[])):
            row=f._paper_search_case(.1,.6,config,(1e-5,),True,1.,.02,2.)
        self.assertFalse(row.nonexistence_supported)
        self.assertFalse(row.diagnostic_complete)
        self.assertEqual(row.status,'ambiguous')

    def test_complete_scan_cannot_override_incomplete_beyn_diagnostic(self):
        config=replace(f.CONFIG,run_physical_diagnostics=False)
        def discover(e,c,margin,nq):
            diag=diagnostic(estimated_rank=0,threshold_rank=0,quadrature_points=nq)
            diag._s0=np.full((2,2),nq/1000.,complex);diag._s1=diag._s0
            return c,np.array([]),diag,[],[]
        scan=f.IndependentScanSummary(.1,.6,'test',1e-6,1.,160,24,0,0,0,0,math.nan,math.nan,False,True,'no-resolved-interior-mode',True)
        with patch.object(f,'_paper_discover_at_margin',discover),patch.object(
            f,'independent_whole_band_svd_scan',return_value=(scan,[])
        ):
            row=f._paper_search_case(.1,.6,config,(1e-5,),True,1.,.02,2.)
        self.assertFalse(row.nonexistence_supported)
        self.assertEqual(row.status,'ambiguous')

    def test_signed_cache_roundtrip_and_configuration_change(self):
        config=replace(f.CONFIG,paper_a_values=(.6,),paper_epsilon_values=(.1,))
        row=f._paper_empty_row(.1,.6,config,'ambiguous','regression')
        with tempfile.TemporaryDirectory() as folder:
            directory=Path(folder)
            with patch.object(f,'run_paper_a_sweep',return_value=[row]) as compute:
                f.load_or_run_paper_a_sweep(config,directory)
            compute.assert_called_once()
            with patch.object(f,'run_paper_a_sweep',side_effect=AssertionError('cache should be valid')):
                loaded=f.load_or_run_paper_a_sweep(config,directory)
                self.assertEqual(loaded[0].assembly_version,f.ASSEMBLY_VERSION)
            with patch.object(f,'run_paper_a_sweep',return_value=[row]) as compute:
                f.load_or_run_paper_a_sweep(replace(config,harmonic_order=40),directory)
            compute.assert_called_once()
            csv_path=directory/'paper_a_sweep_M24.csv'
            csv_path.write_text(csv_path.read_text()+'\n')
            with patch.object(f,'run_paper_a_sweep',return_value=[row]) as compute:
                f.load_or_run_paper_a_sweep(replace(config,harmonic_order=40),directory)
            compute.assert_called_once()


if __name__=='__main__':
    unittest.main()
