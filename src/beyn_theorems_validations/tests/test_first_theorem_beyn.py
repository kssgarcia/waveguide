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


class RunVersioningTests(unittest.TestCase):
    def test_only_pending_selects_stages_without_running_baseline(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg=replace(f.CONFIG,output_directory=tmp)
            with patch.object(f,'CONFIG',cfg),patch.object(f,'main') as execute:
                f.command_line(['--only-pending'])
            active=execute.call_args.args[0]
            self.assertFalse(active.run_main_epsilon_validation)
            self.assertFalse(active.run_internal_convergence_study)
            self.assertTrue(active.run_paper_a_sweep)
            self.assertTrue(active.run_critical_height_study)
            self.assertIn('_pending_',Path(active.output_directory).name)
            scope=json.loads((Path(active.output_directory)/'execution_scope.json').read_text())
            self.assertFalse(scope['prior_results_imported'])

    def test_default_command_keeps_full_execution(self):
        with patch.object(f,'main') as execute:
            f.command_line([])
        execute.assert_called_once_with()

    def test_each_invocation_is_unique_and_preserves_previous_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=replace(f.CONFIG,output_directory=tmp)
            first=f.new_run_config(base)
            sentinel=Path(first.output_directory)/'summary.csv'
            sentinel.write_text('original result')
            second=f.new_run_config(base)
            self.assertNotEqual(first.output_directory,second.output_directory)
            self.assertEqual(sentinel.read_text(),'original result')
            pointer=json.loads((Path(tmp)/'latest_run.json').read_text())
            self.assertEqual(pointer['directory'],second.output_directory)
            manifest=json.loads((Path(first.output_directory)/'run_manifest.json').read_text())
            self.assertEqual(manifest['stage'],'created')

    def test_gate_configuration_runs_in_same_directory_and_failure_is_marked(self):
        with tempfile.TemporaryDirectory() as tmp:
            active=f.new_run_config(replace(f.CONFIG,output_directory=tmp),'gated')
            with patch.object(f,'_run_full',side_effect=RuntimeError('failed solve')):
                with self.assertRaises(RuntimeError):f.main(active)
            self.assertEqual(len([p for p in Path(tmp).iterdir() if p.is_dir()]),1)
            manifest=json.loads((Path(active.output_directory)/'run_manifest.json').read_text())
            self.assertEqual(manifest['stage'],'failed')


class PerformanceSafetyTests(unittest.TestCase):
    def test_lattice_matches_frozen_reference(self):
        import importlib.util
        path = Path(__file__).parent / 'reference/lattice_sums_before_optimization.py'
        spec = importlib.util.spec_from_file_location('reference_lattice', path)
        reference = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reference)
        for k, beta, terms, order in ((.8, 0., 25, 5),(.8+.03j,0.,200,20),
                                      (1.570796-1e-7+.00001j,0.,200,40),(.9,.1,40,8)):
            expected = reference.lattice_sums(4.,k,beta,terms,order)
            actual = f.lattice.lattice_sums(4.,k,beta,terms,order)
            np.testing.assert_array_equal(actual, expected)

    def test_geometry_cache_key_and_immutable_arrays(self):
        from dataclasses import replace
        f.assembly_geometry.cache_clear()
        f.assemble_matrix(.8+.02j,.1,24,f.CONFIG)
        info = f.assembly_geometry.cache_info()
        f.assemble_matrix(.9+.03j,.1,24,f.CONFIG)
        self.assertEqual(f.assembly_geometry.cache_info().hits, info.hits+1)
        for changes in ({'a':.55},{'b':1.1},{'harmonic_order':40},{'finite_difference_step':2e-6}):
            f.assemble_matrix(.8+.02j,.1,24,replace(f.CONFIG,**changes))
        self.assertEqual(f.assembly_geometry.cache_info().misses,5)
        arrays=f.assembly_geometry(.1,24,1.,.6,1e-6,20)
        with self.assertRaises(ValueError): arrays[0][0,0,0]=1.

    def test_checkpoint_reuse_tamper_configuration_and_failure(self):
        from dataclasses import fields
        values={}
        for field in fields(f.ValidationResult):
            kind=str(field.type)
            values[field.name] = ('value' if kind=='str' else False if kind.startswith('bool') else 1 if kind=='int' else .1)
        values.update(epsilon=.1,assembly_version=f.ASSEMBLY_VERSION)
        result=(f.ValidationResult(**values),[diagnostic()],[],[],[],[],np.array([1.,1e-12]),[],float('nan'),None)
        with tempfile.TemporaryDirectory() as tmp:
            cfg=replace(f.CONFIG,output_directory=tmp)
            with patch.object(f,'validate_epsilon',return_value=result) as calculate:
                first=f.baseline_case(.1,cfg);second=f.baseline_case(.1,cfg)
                self.assertEqual(calculate.call_count,1)
                self.assertEqual(second[0],first[0])
                np.testing.assert_array_equal(second[6],first[6])
                path=Path(tmp)/'baseline_cases/epsilon_0.10000000/completed.json'
                envelope=json.loads(path.read_text());envelope['payload']+=' '
                path.write_text(json.dumps(envelope))
                f.baseline_case(.1,cfg);self.assertEqual(calculate.call_count,2)
                f.baseline_case(.1,replace(cfg,harmonic_order=40));self.assertEqual(calculate.call_count,3)
            path.unlink()
            with patch.object(f,'validate_epsilon',side_effect=RuntimeError('refinement failure')):
                with self.assertRaises(RuntimeError): f.baseline_case(.1,cfg)
            self.assertFalse(path.exists())


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


class DeterminantExclusionTests(unittest.TestCase):
    def test_constant_invertible_matrix_has_zero_count(self):
        with patch.object(f,'assemble_matrix',return_value=np.eye(24,dtype=complex)):
            result=f.determinant_contour_exclusion(.02,f.CONFIG)
        self.assertTrue(result['supported'])
        self.assertEqual(len(result['history']),2)

    def test_one_zero_never_supports_absence(self):
        _,_,center,_,_=f.contour_geometry(f.CONFIG)
        def matrix(z,*args):
            A=np.eye(24,dtype=complex);A[0,0]=z-center;return A
        with patch.object(f,'assemble_matrix',side_effect=matrix):
            result=f.determinant_contour_exclusion(.02,f.CONFIG)
        self.assertTrue(result['complete'])
        self.assertFalse(result['supported'])
        self.assertEqual(result['history'][-1]['count'],1)

    def test_nan_or_failure_never_supports_exclusion(self):
        for value in (np.full((24,24),np.nan),np.zeros((24,24))):
            with patch.object(f,'assemble_matrix',return_value=value):
                result=f.determinant_contour_exclusion(.02,f.CONFIG)
            self.assertFalse(result['supported'])
            self.assertFalse(result['complete'])
        with patch.object(f,'assemble_matrix',side_effect=ValueError('injected')):
            self.assertFalse(f.determinant_contour_exclusion(.02,f.CONFIG)['supported'])

    def test_ill_conditioned_matrix_never_supports_exclusion(self):
        A=np.eye(24,dtype=complex);A[0]=A[1];A[0,0]=1e-12
        with patch.object(f,'assemble_matrix',return_value=A):
            result=f.determinant_contour_exclusion(.02,f.CONFIG)
        self.assertFalse(result['supported'])

    def test_phase_underresolution_never_supports_absence(self):
        cfg=replace(f.CONFIG,nonexistence_max_phase_step=-1.)
        with patch.object(f,'assemble_matrix',return_value=np.eye(24,dtype=complex)):
            self.assertFalse(f.determinant_contour_exclusion(.02,cfg)['supported'])

class ScanEndpointTests(unittest.TestCase):
    def test_low_k_endpoint_is_evaluated_at_exact_configured_k(self):
        cfg=f.CONFIG
        endpoint=f.sigma_from_kb(cfg.beyn_low_k_margin,cfg)
        with patch.object(f,'singular_metrics',return_value=(.5,1.,.5)) as metrics:
            self.assertEqual(f.absolute_singular_value_from_sigma(endpoint,.02,24,cfg),.5)
            self.assertEqual(metrics.call_args.args[0],cfg.beyn_low_k_margin)
            self.assertEqual(f.relative_singular_value_from_sigma(endpoint,.02,24,cfg),.5)

    def test_cutoff_and_outside_low_k_remain_excluded(self):
        cfg=f.CONFIG
        self.assertTrue(math.isinf(f.absolute_singular_value_from_sigma(0.,.02,24,cfg)))
        outside=f.sigma_from_kb(.5*cfg.beyn_low_k_margin,cfg)
        self.assertTrue(math.isinf(f.absolute_singular_value_from_sigma(outside,.02,24,cfg)))

    def test_endpoint_failure_is_not_hidden(self):
        endpoint=f.sigma_from_kb(f.CONFIG.beyn_low_k_margin,f.CONFIG)
        with patch.object(f,'singular_metrics',side_effect=ValueError('failure')):
            with self.assertRaises(ValueError):
                f.absolute_singular_value_from_sigma(endpoint,.02,24,f.CONFIG)

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
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        patcher = patch.object(f, 'CONFIG', replace(f.CONFIG, output_directory=temporary.name))
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_new_audit_fields_survive_csv_roundtrip(self):
        cfg=f.CONFIG
        row=f._paper_empty_row(.1,.6,cfg,'ambiguous','test')
        row.beyn_threshold_rank=2;row.beyn_gap_rank=1
        row.beyn_s0_relative_change=.002;row.beyn_s1_relative_change=.003
        row.physical_checks_run=True;row.physical_checks_passed=False
        row.wall_relative_residual=.01
        path=Path(cfg.output_directory)/'rows.csv'
        f.write_dataclass_csv(path,[row])
        saved=f.read_paper_sweep_csv(path)[0]
        self.assertEqual(saved.beyn_threshold_rank,2)
        self.assertEqual(saved.beyn_gap_rank,1)
        self.assertEqual(saved.beyn_s1_relative_change,.003)
        self.assertTrue(saved.physical_checks_run)
        self.assertFalse(saved.physical_checks_passed)
        self.assertEqual(saved.wall_relative_residual,.01)

    def test_evaluated_rejected_valley_requires_independent_exclusion(self):
        def discover(e,c,margin,nq):
            diag=diagnostic(estimated_rank=0,threshold_rank=0,gap_rank=0,quadrature_points=nq)
            diag._s0=np.zeros((2,2),complex);diag._s1=diag._s0
            return c,np.array([]),diag,[],[]
        mode=f.ModeResult(.1,1,'scan',1.,0.,1.,.02,.1,2.,.05,1.,math.nan,False,.05,.05,.1)
        scan=f.IndependentScanSummary(.1,.35,'test',1e-6,1.,160,24,1,1,0,1,.05,1.,False,False,'ambiguous',False)
        for independent_zero in (True,False):
            with patch.object(f,'_paper_discover_at_margin',discover),patch.object(
                f,'independent_whole_band_svd_scan',return_value=(scan,[mode])
            ),patch.object(f,'determinant_contour_exclusion',return_value={'supported':independent_zero}):
                row=f._paper_search_case(.1,.35,replace(f.CONFIG,nonexistence_refinement_factors=()),(1e-5,),False,math.nan,math.nan,-1.)
            self.assertEqual(row.nonexistence_supported,independent_zero)

    def test_nan_candidate_cannot_be_closed_by_exclusion(self):
        mode=f.ModeResult(.1,1,'scan',1.,0.,1.,.02,.1,2.,.05,1.,math.nan,False)
        self.assertFalse(f.mode_evaluation_complete(mode))

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


class ExtendedFieldTests(unittest.TestCase):
    def test_modal_field_outside_lattice_disk_and_near_cutoff(self):
        for a, eps in ((.5, .449474), (.5, .474737), (.55, .449474)):
            cfg=replace(f.CONFIG,a=a)
            for kb in (1., np.pi/2-1e-7):
                green=f.make_extended_field_green(kb,cfg)
                for x in (-3.,3.):
                    for y in (-cfg.b-a, -.2-a, cfg.b-a):
                        actual=green(x,y,-eps,eps)
                        expected=modal_green(x+eps,y+cfg.b+a,eps+cfg.b+a,kb)
                        self.assertTrue(np.isfinite(actual))
                        self.assertLess(abs(actual-expected),1e-12*max(1.,abs(expected)))

    def test_modal_and_lattice_agree_in_overlap(self):
        cfg=replace(f.CONFIG,harmonic_order=40)
        lattice_green,_=f.make_green_functions(1.,cfg)
        extended=f.make_extended_field_green(1.,cfg)
        for x in (-1.,1.,1.2):
            actual=extended(x,-.3,0.,.1)
            expected=lattice_green(x,-.3,0.,.1)
            self.assertLess(abs(actual-expected),1e-11)

if __name__=='__main__':
    unittest.main()
