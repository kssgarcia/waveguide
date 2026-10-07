"""Focused convergence evidence, gated by the current representative run.

Runs the remaining baseline epsilons and sampled subcritical geometries.
Does not rerun the 99-point publication sweep; never reuses unsigned evidence.
"""
from __future__ import annotations
import contextlib
import csv
import json
import multiprocessing as mp
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import first_theorem_beyn as f
from validate_first_theorem_beyn import representative_case, legacy_snapshot, run_gate


def subcritical_case(epsilon: float, config: f.Config) -> dict:
    root = Path(config.output_directory) / 'convergence_study'
    with (root / f'subcritical_{epsilon:.2f}.log').open('w', buffering=1) as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        row = f.evaluate_paper_sweep_geometry(epsilon, .35, config)
    return asdict(row)


def main() -> None:
    config = f.new_run_config(f.CONFIG, "convergence")
    if run_gate(config):
        raise RuntimeError("Representative gate failed; convergence study not started")
    root = Path(config.output_directory)
    gate = json.loads((root / 'validation_gate.json').read_text())
    if not gate['passed'] or gate['configuration_fingerprint'] != f.configuration_fingerprint(config):
        raise RuntimeError('Run the current representative gate before this study')
    before = legacy_snapshot()
    output = root / 'convergence_study'
    output.mkdir(exist_ok=True)
    f.write_run_manifest(output, config, 'convergence-study-running')
    results = list(gate['cases'])
    subcritical = []
    remaining = [e for e in config.epsilon_values if e not in (.02, .10, .20)]
    with ProcessPoolExecutor(max_workers=4, mp_context=mp.get_context('spawn')) as pool:
        jobs = {pool.submit(representative_case, e, config): ('baseline', e) for e in remaining}
        jobs.update({pool.submit(subcritical_case, e, config): ('subcritical', e) for e in (.02, .04, .06, .10)})
        for job in as_completed(jobs):
            kind, epsilon = jobs[job]
            row = job.result()  # Exceptions stop completion; no success manifest.
            (results if kind == 'baseline' else subcritical).append(row)
            print(json.dumps({'kind': kind, **row}), flush=True)
    results.sort(key=lambda r: r['epsilon'])
    subcritical.sort(key=lambda r: r['epsilon'])
    summaries = []
    for e in config.epsilon_values:
        path = root / 'representative' / f'epsilon_{e:.2f}' / 'summary.csv'
        summaries.extend(csv.DictReader(path.open()))
    with (output / 'baseline_summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader(); writer.writerows(summaries)
    assert before == legacy_snapshot(), 'Legacy changed'
    report = {'configuration_fingerprint': f.configuration_fingerprint(config),
              'baseline': results, 'subcritical': subcritical,
              'legacy_unchanged': True, 'M': config.bem_M,
              'scope': '10 baseline and 4 subcritical cases; no 99-point sweep or critical-height rerun'}
    (output / 'results.json').write_text(json.dumps(report, indent=2)+'\n')
    f.write_run_manifest(output, config, 'convergence-study-completed')
    f.write_run_manifest(root, config, 'convergence-study-completed')


if __name__ == '__main__':
    main()
