"""Gated v9 validation. No long sweep starts until all regressions/case gates pass.

Usage: python validate_first_theorem_beyn.py [--full-after-gate]
The flag explicitly authorizes the production sweep after the representative gate.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import multiprocessing as mp
import sys
import unittest
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import first_theorem_beyn as f


def legacy_snapshot() -> dict[str, str]:
    return {
        str(path.relative_to(f.LEGACY_OUTPUT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(f.LEGACY_OUTPUT.rglob('*')) if path.is_file()
    }


def representative_case(epsilon: float, config: f.Config) -> dict:
    directory=Path(config.output_directory)/'representative'/f'epsilon_{epsilon:.2f}'
    f.write_run_manifest(directory,config,'representative-running')
    with (directory/'run.log').open('w',buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        result,diags,eigen,convergence,refinements,modes,_,_,_,physical=f.validate_epsilon(epsilon,config)
        f.write_dataclass_csv(directory/'summary.csv',[result])
        f.write_dataclass_csv(directory/'beyn_diagnostics.csv',diags)
        f.write_dataclass_csv(directory/'beyn_convergence.csv',convergence)
        f.write_dataclass_csv(directory/'beyn_raw_eigenvalues.csv',eigen)
        f.write_dataclass_csv(directory/'mode_refinement.csv',refinements)
        f.write_dataclass_csv(directory/'mode_results.csv',modes)
        f.write_dataclass_csv(directory/'physical_diagnostics.csv',[physical] if physical else [])
        f.print_result(result,config)
        shift=float('nan')
        high_mode=None
        if result.mode_exists_supported:
            high_config=replace(config,harmonic_order=40)
            seed=f.DiscoverySeed(result.kb_numerical,0.,'representative-H40')
            _,high_mode=f.run_candidate_refinement(
                epsilon,1,seed,(.8*result.sigma_numerical,1.2*result.sigma_numerical),high_config
            )
            shift=abs(high_mode.sigma_numerical-result.sigma_numerical)/abs(high_mode.sigma_numerical)
            f.write_dataclass_csv(directory/'harmonic_order_40_mode.csv',[high_mode])
        passed=bool(
            config.bem_M==24 and result.mode_exists_supported
            and high_mode is not None and high_mode.resolved
            and np.isfinite(shift) and shift<=config.harmonic_stability_tolerance
            and (not result.uniqueness_supported or (
                result.beyn_moments_converged and result.beyn_rank_stable
                and result.tighter_margin_rank_consistent is True
            ))
        )
        report={
            'epsilon':epsilon,'passed':passed,'mode_exists_supported':result.mode_exists_supported,
            'uniqueness_supported':result.uniqueness_supported,
            'uniqueness_status':result.uniqueness_status,'sigma':result.sigma_numerical,
            'harmonic_20_vs_40_relative_sigma_shift':shift,
            'normalized_residual':result.normalized_residual,
            'equilibrated_residual':result.equilibrated_residual,
            'beyn_s0_change':result.beyn_final_s0_relative_change,
            'beyn_s1_change':result.beyn_final_s1_relative_change,
        }
        (directory/'gate.json').write_text(json.dumps(report,indent=2)+'\n')
        f.write_run_manifest(directory,config,'representative-passed' if passed else 'representative-failed')
        return report


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full-after-gate',action='store_true')
    args=parser.parse_args()
    config=f.CONFIG
    if config.bem_M != 24:
        raise ValueError('This corrected validation is fixed at M=24')
    before=legacy_snapshot()
    directory=Path(config.output_directory)
    f.write_run_manifest(directory,config,'regressions-running')
    suite=unittest.defaultTestLoader.discover(str(Path(__file__).resolve().parents[1] / 'tests'),pattern='test_first_theorem_beyn.py')
    with (directory/'regression_tests.log').open('w',buffering=1) as stream:
        tests=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    if not tests.wasSuccessful():
        f.write_run_manifest(directory,config,'regressions-failed-no-full-run')
        return 1
    print(f'Regressions passed: {tests.testsRun}',flush=True)
    reports=[]
    with ProcessPoolExecutor(max_workers=min(3,f.paper_worker_count(config)),mp_context=mp.get_context('spawn')) as pool:
        futures={pool.submit(representative_case,e,config):e for e in (.02,.10,.20)}
        for future in as_completed(futures):
            e=futures[future]
            try:
                report=future.result()
            except Exception as exc:
                report={'epsilon':e,'passed':False,'error':f'{type(exc).__name__}: {exc}'}
            reports.append(report)
            print(json.dumps(report),flush=True)
    reports.sort(key=lambda row:row['epsilon'])
    passed=all(row['passed'] for row in reports) and before==legacy_snapshot()
    gate={
        'passed':passed,'regression_tests':tests.testsRun,'cases':reports,
        'legacy_unchanged':before==legacy_snapshot(),
        'assembly_version':f.ASSEMBLY_VERSION,'configuration_fingerprint':f.configuration_fingerprint(config),
        'gate_definition':'All regressions; local existence + physical checks at M24; H20/H40 stable. Unconverged uniqueness remains ambiguous, not a failed existence claim.',
    }
    (directory/'validation_gate.json').write_text(json.dumps(gate,indent=2)+'\n')
    (directory/'legacy_sha256.json').write_text(json.dumps(before,indent=2)+'\n')
    if not passed:
        f.write_run_manifest(directory,config,'representative-gate-failed-no-full-run')
        print('Gate failed: production sweep was NOT started.',flush=True)
        return 1
    f.write_run_manifest(directory,config,'representative-gate-passed')
    if args.full_after_gate:
        print('Gate passed; starting full corrected validation.',flush=True)
        try:
            f.main()
        except Exception:
            f.write_run_manifest(directory,config,'full-validation-failed')
            raise
        if before != legacy_snapshot():
            raise RuntimeError('Legacy evidence changed during corrected validation')
    return 0


if __name__=='__main__':
    sys.exit(main())
