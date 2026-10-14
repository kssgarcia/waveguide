"""Resume an existing compatible run, preserving completed scientific work.

Usage: python runners/resume_validation.py --theorem 2.1 --run PATH
No changes to assembly, search, tolerances or original fingerprint definitions.
"""
import argparse,csv,json,sys,hashlib,os,fcntl
from pathlib import Path
from datetime import datetime
from contextlib import redirect_stdout,redirect_stderr
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import first_theorem_beyn as f
import second_theorem_X as x
import second_theorem_Y as y
import execution_harness as h
from refresh_reviewed_figures import records


def load_config(module,run):
    run=Path(run).resolve();manifest=json.loads((run/'run_manifest.json').read_text())
    cfg=module.Config(**{k:tuple(v) if isinstance(v,list) else v for k,v in manifest['config'].items()})
    if Path(cfg.output_directory).resolve()!=run:raise ValueError('Run path differs from saved configuration')
    expected=f.configuration_fingerprint(cfg) if module is f else h.fingerprint(module,cfg)
    saved=manifest.get('configuration_fingerprint') if module is f else manifest.get('fingerprint')
    if expected!=saved:raise ValueError('Source/configuration changed: cannot reuse this run; start a new execution')
    return cfg,manifest


def valid_paper(row,a,epsilon,cfg):
    return (isinstance(row,f.PaperSweepRow) and row.a==a and row.epsilon==epsilon
            and row.M==cfg.bem_M and row.assembly_version==f.ASSEMBLY_VERSION
            and row.status!='error')


def paper_checkpoint_worker(payload):
    a,epsilon,cfg=payload
    directory=Path(cfg.output_directory)/'paper_cases'/f'a_{a:.12g}_epsilon_{epsilon:.12g}'
    return h.checkpoint(directory,f.configuration_fingerprint(cfg,epsilon,a),True,
        lambda:f.evaluate_paper_sweep_geometry(epsilon,a,cfg),h.encode,lambda p:h.decode(p,f),
        lambda r:valid_paper(r,a,epsilon,cfg))


def prepare_first(cfg):
    run=Path(cfg.output_directory);adopted=[];missing=[]
    # Older runs have atomically written geometry certificates before a final
    # CSV. Adopt these ONLY after the outer run fingerprint was verified.
    for a in cfg.paper_a_values:
        for epsilon in f.admissible_paper_epsilon_values(a,cfg):
            key=f'a_{a:.12g}_epsilon_{epsilon:.12g}'
            checkpoint=run/'paper_cases'/key/'completed.json'
            source=run/'geometry_diagnostics'/key/'certificate.json'
            if checkpoint.exists():continue
            try:
                row=f.PaperSweepRow(**json.loads(source.read_text()))
                if not valid_paper(row,a,epsilon,cfg):raise ValueError('Invalid geometry certificate')
                h.checkpoint(checkpoint.parent,f.configuration_fingerprint(cfg,epsilon,a),True,
                    lambda:row,h.encode,lambda p:h.decode(p,f),lambda r:valid_paper(r,a,epsilon,cfg))
                adopted.append({'key':key,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()})
            except (OSError,ValueError,TypeError):missing.append(key)
    h.atomic_json(run/'resume_recovery.json',{'outer_fingerprint_verified':True,'adopted_cases':adopted,'pending_cases':missing})
    print(f'Recovered {len(adopted)} geometry certificates; pending {len(missing)}',flush=True)
    f._paper_worker=paper_checkpoint_worker
    original_internal=f.run_internal_convergence_study
    def internal(results,c):
        path=run/'internal_convergence.csv'
        if path.exists():
            try:rows=records(path,f.InternalConvergenceRow)
            except (OSError,ValueError,TypeError,KeyError):rows=[]
            expected={(name,float(v)) for name,values in [('finite_difference_step',c.finite_difference_steps_test),('lattice_terms',c.lattice_terms_test),('harmonic_order',c.harmonic_orders_test)] for v in values}
            if {(r.parameter,r.value) for r in rows}==expected:
                print('Reusing completed internal study (run fingerprint verified)',flush=True);return rows
        return original_internal(results,c)
    original_critical=f.run_critical_height_study
    def critical(c):
        if (run/'critical_height.csv').exists() and (run/'critical_height_probes.csv').exists():
            try:
                rows=records(run/'critical_height.csv',f.CriticalHeightRow);probes=records(run/'critical_height_probes.csv',f.CriticalHeightProbeRow)
            except (OSError,ValueError,TypeError,KeyError):rows=[];probes=[]
            if {r.epsilon for r in rows}==set(c.critical_height_epsilon_values) and {p.epsilon for p in probes}==set(c.critical_height_epsilon_values) and all(r.status=='bracketed' for r in rows):
                print('Reusing completed critical-height study (run fingerprint verified)',flush=True);return rows,probes
        return original_critical(c)
    f.run_internal_convergence_study=internal;f.run_critical_height_study=critical


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--theorem',choices=('2.1','X','Y'),required=True)
    parser.add_argument('--run',type=Path,default=Path('latest'),help='Saved run path; default: latest run of this theorem')
    args=parser.parse_args();module={'2.1':f,'X':x,'Y':y}[args.theorem]
    selected=args.run
    if selected==Path('latest'):
        selected=Path(json.loads((Path(module.CONFIG.output_directory)/'latest_run.json').read_text())['directory'])
    cfg,prior=load_config(module,selected);run=Path(cfg.output_directory)
    # A lock prevents concurrent writers. On interruption it is released.
    lock=run/'.resume.lock'
    lock_handle=lock.open('a+')
    try:fcntl.flock(lock_handle.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        lock_handle.close();raise RuntimeError('Another resume process owns this run')
    lock_handle.seek(0);lock_handle.truncate();lock_handle.write(str(os.getpid()));lock_handle.flush()
    event=run/'resume_history'/datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    event.mkdir(parents=True,exist_ok=False)
    h.atomic_json(event/'original_manifest.json',prior)
    h.atomic_json(event/'resume.json',{'theorem':args.theorem,'run':str(run),'status':'running','python':sys.version})
    try:
        with (event/'execution.log').open('w',buffering=1) as log,redirect_stdout(log),redirect_stderr(log):
            print(f'Resuming {args.theorem}, unchanged fingerprint; original files preserved in their run',flush=True)
            if module is f:
                prepare_first(cfg);f.main(cfg)
            elif module is x:
                from validate_second_theorem_X import run as operation
                h.manifest(module,cfg,'running');operation(cfg);h.manifest(module,cfg,'completed')
            else:
                from validate_second_theorem_Y import check_results,check_internal_sensitivity,write_verification_report
                h.manifest(module,cfg,'running');module._run_validation(cfg)
                passed=check_results(run) and check_internal_sensitivity(run)
                summaries=list(csv.DictReader((run/'summary.csv').open()))
                all_modes=bool(len(summaries)==len(cfg.epsilon_values)
                               and all(r['mode_exists_supported']=='True' for r in summaries))
                h.atomic_json(run/'verification_gate.json',{'representative_independent_checks_passed':passed,
                    'all_requested_modes_supported':all_modes,'resumed':True})
                write_verification_report(run)
                if not passed or not all_modes:raise RuntimeError('Independent verification failed')
                h.manifest(module,cfg,'completed')
        h.atomic_json(event/'resume.json',{'status':'completed','run':str(run),'theorem':args.theorem})
        print(f'Completed: {run}',flush=True)
    except BaseException as exc:
        state='interrupted' if isinstance(exc,KeyboardInterrupt) else 'failed'
        (f.write_run_manifest(run,cfg,state) if module is f else h.manifest(module,cfg,state))
        h.atomic_json(event/'resume.json',{'status':state,'error':f'{type(exc).__name__}: {exc}'})
        raise
    finally:
        fcntl.flock(lock_handle.fileno(),fcntl.LOCK_UN);lock_handle.close()

if __name__=='__main__':main()
