"""Recompute selected saved-run statuses, with fresh versioned output."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import first_theorem_beyn as f


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_run', type=Path)
    parser.add_argument('--status', choices=('error', 'ambiguous'), default='error')
    parser.add_argument('--critical-height', action='store_true')
    args = parser.parse_args()
    source = args.source_run.resolve()
    csv_path = source / 'paper_a_sweep_M24.csv'
    rows = [r for r in f.read_paper_sweep_csv(csv_path) if r.status == args.status]
    if not rows:
        raise SystemExit(f'No {args.status} rows to retry')
    saved = json.loads((source / 'run_manifest.json').read_text())['config']
    # JSON arrays are restored to the immutable tuple configuration fields.
    saved = {k: tuple(v) if isinstance(v, list) else v for k, v in saved.items()}
    cfg = f.Config(**saved)
    if args.status == 'ambiguous':
        cfg = replace(cfg, nonexistence_scan_max_candidates=max(
            cfg.nonexistence_scan_max_candidates, f.CONFIG.nonexistence_scan_max_candidates))
    cfg = replace(cfg, output_directory=f.CONFIG.output_directory)
    active = f.new_run_config(cfg, 'resolve-pending' if args.status == 'ambiguous' else 'retry-errors')
    output = Path(active.output_directory)
    f.harness.atomic_json(output / 'execution_scope.json', {
        'source_run': str(source),
        'source_csv_sha256': hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        'scope': f'only prior {args.status} rows; fresh numerical evaluations',
        'critical_height': args.critical_height,
        'candidate_budget': active.nonexistence_scan_max_candidates,
        'source_candidate_budget': saved['nonexistence_scan_max_candidates'],
        'cases': [{'a': r.a, 'epsilon': r.epsilon} for r in rows],
        'prior_results_imported': False,
    })
    f.write_run_manifest(output, active, 'targeted-validation-running')
    print(f'Output: {output}', flush=True)
    completed = []
    try:
        with ProcessPoolExecutor(max_workers=min(3, len(rows))) as pool:
            futures = {pool.submit(f._paper_worker, (r.a, r.epsilon, active)): r for r in rows}
            for future in as_completed(futures):
                original = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    result = f._paper_empty_row(original.epsilon, original.a, active,
                                               'error', f'{type(exc).__name__}: {exc}')
                completed.append(result)
                f.write_dataclass_csv(output / 'paper_a_sweep_M24.csv',
                                      sorted(completed, key=lambda r: (r.a, r.epsilon)))
                print(f'a={result.a} epsilon={result.epsilon}: {result.status}; '
                      f'physical={result.physical_checks_passed}; '
                      f'unique={result.uniqueness_supported}', flush=True)
        if args.critical_height:
            critical, probes = f.run_critical_height_study(active)
            f.write_dataclass_csv(output / 'critical_height.csv', critical)
            f.write_dataclass_csv(output / 'critical_height_probes.csv', probes)
        failures = sum(r.status == 'error' for r in completed)
        f.write_run_manifest(output, active, 'targeted-validation-completed' if not failures else 'targeted-validation-with-errors')
        lines = ['# Targeted Theorem 2.1 validation', '',
                 f'Source run: `{source.name}`. Original files preserved.',
                 f'Assembly/Green version: `{f.ASSEMBLY_VERSION}`; M={active.bem_M}.',
                 f'Only prior {args.status} geometries were recomputed. Critical study: {args.critical_height}. No full sweep.', '',
                 '| a | epsilon | status | mode supported | uniqueness supported | physical passed |',
                 '|---|---|---|---|---|---|']
        lines += [f'| {r.a} | {r.epsilon} | {r.status} | {r.mode_exists_supported} | '
                  f'{r.uniqueness_supported} | {r.physical_checks_passed} |' for r in completed]
        (output / 'results_audit.md').write_text('\n'.join(lines)+'\n')
        if failures:
            raise SystemExit(1)
    except BaseException:
        f.write_run_manifest(output, active, 'targeted-validation-interrupted')
        raise


if __name__ == '__main__':
    main()
