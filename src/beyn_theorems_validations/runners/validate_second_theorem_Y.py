"""Versioned Y smoke/full run plus independent modal/operator and field checks."""
import sys,argparse,csv,json
import hashlib
from pathlib import Path
from dataclasses import replace,asdict
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import second_theorem_Y as y
from check_X_modal_operator import modal_matrix


def check_results(run):
    run=Path(run)
    manifest=json.loads((run/'run_manifest.json').read_text())
    cfg=y.Config(**{k:tuple(v) if isinstance(v,list) else v for k,v in manifest['config'].items()})
    if manifest['fingerprint']!=y.harness.fingerprint(y,cfg):
        raise RuntimeError('Source/configuration changed after execution; independent check cannot silently certify it')
    diagnostics=[]
    for epsilon in (.02,.10,.20):
        directory=run/'certificates'/f'epsilon_{epsilon:.8f}'
        if not directory.exists():continue
        mode=y.TunedMode(**json.loads((directory/'mode.json').read_text()))
        k=mode.kb_numerical;a=mode.a_numerical
        A=y.assemble_full_matrix(k,epsilon,a,24,cfg)
        B=modal_matrix(k,epsilon,cfg,12000,geometry=lambda t,e,c:y.obstacle_geometry(t,e,a,c))
        C=modal_matrix(k,epsilon,cfg,24000,geometry=lambda t,e,c:y.obstacle_geometry(t,e,a,c))
        _,_,_,_,v=y.full_singular_pair(k,epsilon,a,24,cfg)
        reference=y.physical_mode_diagnostics(epsilon,mode,cfg)
        variants=[]
        for label,c in [('quadrature-64',replace(cfg,physical_quadrature_points=64)),
                        ('quadrature-128',replace(cfg,physical_quadrature_points=128)),
                        ('source-oversample-8',replace(cfg,physical_boundary_source_oversample_factor=8)),
                        ('sections-shifted',replace(cfg,physical_x_over_b=(1.6,2.1,2.6,3.1)))]:
            p=y.physical_mode_diagnostics(epsilon,mode,c)
            variants.append({'variant':label,'physical_supported':y.physical_is_supported(p),
                'BIE_residual':p.boundary_integral_relative_residual,
                'open_channel_ratio':p.first_open_channel_relative_amplitude,
                'modal_projection_discrepancy':p.modal_projection_discrepancy,
                'decay_relative_error':max(p.decay_relative_error_left,p.decay_relative_error_right)})
        H40=y.assemble_full_matrix(k,epsilon,a,24,replace(cfg,harmonic_order=40))
        result={'epsilon':epsilon,'M':24,'mode_exists_supported':mode.resolved and y.physical_is_supported(reference),
            'modal_terms':[12000,24000],
            'matrix_relative_difference':float(np.linalg.norm(A-C)/np.linalg.norm(C)),
            'modal_truncation_relative_change':float(np.linalg.norm(B-C)/np.linalg.norm(C)),
            'modal_vector_residual':float(np.linalg.norm(C@v)/(np.linalg.norm(C,2)*np.linalg.norm(v))),
            'H20_H40_matrix_relative_change':float(np.linalg.norm(A-H40)/np.linalg.norm(H40)),
            'physical_variants':variants}
        result['independent_checks_passed']=bool(result['mode_exists_supported'] and
            result['matrix_relative_difference']<1e-8 and result['modal_truncation_relative_change']<1e-8
            and result['modal_vector_residual']<1e-8 and result['H20_H40_matrix_relative_change']<1e-8
            and all(v['physical_supported'] for v in variants))
        diagnostics.append(result)
        y.harness.atomic_json(run/'independent_validation.json',diagnostics)
        print(json.dumps(result),flush=True)
    unbiased=check_unbiased_placement(run,cfg)
    return bool(len(diagnostics)==3 and all(r['independent_checks_passed'] for r in diagnostics) and unbiased)


def check_unbiased_placement(run,cfg=None):
    """Recover placement without centering its search on formula (2.7)."""
    run=Path(run)
    if cfg is None:
        manifest=json.loads((run/'run_manifest.json').read_text())
        cfg=y.Config(**{k:tuple(v) if isinstance(v,list) else v for k,v in manifest['config'].items()})
        if manifest['fingerprint']!=y.harness.fingerprint(y,cfg):raise RuntimeError('Fingerprint mismatch')
    rows=[]
    for epsilon in (.02,.10):
        mode=y.TunedMode(**json.loads((run/'certificates'/f'epsilon_{epsilon:.8f}'/'mode.json').read_text()))
        result,source,lam,norm=y._solve_complex_eigenvalue_root(epsilon,0.,mode.sigma_asymptotic,32.,
            [('zero-placement',np.array([0.,0.])),('positive-placement',np.array([4.,0.]))],cfg)
        a,sigma,kb=y._joint_candidate_from_variables(epsilon,0.,mode.sigma_asymptotic,*result.x,cfg)
        candidate=replace(mode,a_numerical=a,kb_numerical=kb,sigma_numerical=sigma)
        physical=y.physical_mode_diagnostics(epsilon,candidate,cfg)
        certificate=y.certify_full_mode(kb,epsilon,a,24,True,mode.drop_factor,cfg)
        row={'epsilon':epsilon,'search_a1':0.,'correction_bound':32.,'a_recovered':a,
            'placement_difference':abs(a-mode.a_numerical),'relative_sigma_difference':abs(sigma-mode.sigma_numerical)/mode.sigma_numerical,
            'root_norm':norm,'physical_supported':y.physical_is_supported(physical),
            'spectral_supported':certificate['supported'],'attempts':result.attempts}
        row['passed']=bool(result.success and row['physical_supported'] and row['spectral_supported']
                           and row['placement_difference']<1e-8 and row['relative_sigma_difference']<1e-6)
        rows.append(row)
        y.harness.atomic_json(run/'independent_placement.json',rows)
        print('unbiased placement '+json.dumps(row),flush=True)
    return bool(len(rows)==2 and all(r['passed'] for r in rows))


def check_internal_sensitivity(run):
    """Retain full signed-channel/physical records for already retuned cases."""
    run=Path(run);manifest=json.loads((run/'run_manifest.json').read_text())
    cfg=y.Config(**{k:tuple(v) if isinstance(v,list) else v for k,v in manifest['config'].items()})
    if manifest['fingerprint']!=y.harness.fingerprint(y,cfg):raise RuntimeError('Fingerprint mismatch')
    if not cfg.run_internal_convergence_study:return True
    path=run/'internal_convergence.csv'
    if not path.exists():return False
    rows=list(csv.DictReader(path.open()));output=[]
    summaries=list(csv.DictReader((run/'summary.csv').open()))
    target=min(summaries,key=lambda r:abs(float(r['epsilon'])-cfg.internal_convergence_epsilon))
    eps=float(target['epsilon'])
    baseline=y.TunedMode(**json.loads((run/'certificates'/f'epsilon_{eps:.8f}'/'mode.json').read_text()))
    for row in rows:
        parameter=row['parameter'];value=float(row['value'])
        variant=replace(cfg,**{parameter:int(value) if parameter in ('harmonic_order','lattice_terms') else value})
        record={'parameter':parameter,'value':value,'original_status':row['status'],
                'epsilon':eps,'M':24,'physical':None,'certificate':None,'verified':False}
        try:
            a=float(row['a_numerical']);sigma=float(row['sigma_bem']);kb=float(row['kb'])
            if not np.all(np.isfinite([a,sigma,kb])):raise ValueError('Incomplete retuned endpoint')
            mode=replace(baseline,a_numerical=a,sigma_numerical=sigma,kb_numerical=kb)
            _,_,rel,scaled=y.full_singular_metrics(kb,eps,a,24,variant)
            drop=y._spectral_drop_factor(eps,a,sigma,max(rel,scaled),variant)
            certificate=y.certify_full_mode(kb,eps,a,24,True,drop,variant)
            physical=y.physical_mode_diagnostics(eps,mode,variant)
            record.update(a_numerical=a,sigma=sigma,kb=kb,physical=asdict(physical),certificate=certificate)
            record['verified']=bool(row['status']=='certified' and certificate['supported'] and y.physical_is_supported(physical))
        except Exception as exc:record['error']=f'{type(exc).__name__}: {exc}'
        record['reference_quality_variant']=bool(
            parameter=='harmonic_order' and value>=cfg.harmonic_order or
            parameter=='lattice_terms' and value>=100 or
            parameter=='finite_difference_step' and value in (1e-7,1e-6))
        record['stable']=bool(record['verified'] and float(row['relative_sigma_shift_from_reference'])<5e-4
                              and float(row['a_shift_from_reference'])<1e-6*cfg.b)
        output.append(record)
        y.harness.atomic_json(run/'internal_physical_certificates.json',output)
    required=[r for r in output if r['reference_quality_variant']]
    passed=bool(len(output)==len(rows) and required and all(r['stable'] for r in required))
    y.harness.atomic_json(run/'internal_sensitivity_gate.json',{'passed':passed,
        'required_variants':len(required),'limitation':'Coarse harmonic orders and extreme FD steps are diagnostics; failed cases remain ambiguous.'})
    print(f'Internal sensitivity gate: {passed}, {len(required)} reference-quality variants',flush=True)
    return passed


def write_verification_report(run):
    run=Path(run);manifest=json.loads((run/'run_manifest.json').read_text())
    rows=list(csv.DictReader((run/'summary.csv').open()));physical=list(csv.DictReader((run/'physical_diagnostics.csv').open()))
    phys={float(r['epsilon']):r for r in physical}
    independent=json.loads((run/'independent_validation.json').read_text())
    placement=json.loads((run/'independent_placement.json').read_text())
    verified=sum(r['mode_exists_supported']=='True' for r in rows)
    small=[r for r in rows if float(r['epsilon'])<=.1+1e-15]
    lines=['# Resultado verificado de Y, teorema 2.3(iv)', '',
        f"Ejecución `{run.name}`; implementación `{manifest['assembly_version']}`; M=24.", '',
        f'**Existencia física apoyada numéricamente: {verified}/{len(rows)} casos.**',
        '**Unicidad de la banda completa: ambiguous/not-established.** El screen finito no es un conteo certificado.', '',
        f"Screen refinado aprobado en {sum(r['uniqueness_screen_passed']=='True' for r in small)}/{len(small)} casos pequeños; no significa unicidad global.", '',
        '| epsilon | a | sigma | Residual normalizado | Canal abierto | BIE off-grid | Modo apoyado |',
        '| --- | --- | --- | --- | --- | --- | --- |']
    for row in rows:
        p=phys.get(float(row['epsilon']),{})
        lines.append(f"| {row['epsilon']} | {float(row['a_numerical']):.10g} | {float(row['sigma_numerical']):.10g} | {float(row['normalized_residual']):.3e} | {float(p.get('first_open_channel_relative_amplitude','nan')):.3e} | {float(p.get('boundary_integral_relative_residual','nan')):.3e} | {row['mode_exists_supported']} |")
    lines += ['', '## Contrastes independientes', '',
        f"Matrices modales independientes 12000/24000: diferencia lattice/modal máxima {max(r['matrix_relative_difference'] for r in independent):.3e}; residual del vector en matriz modal máximo {max(r['modal_vector_residual'] for r in independent):.3e}.",
        f"Harmonic order 20/40: cambio relativo máximo {max(r['H20_H40_matrix_relative_change'] for r in independent):.3e}.",
        'Los tres representantes .02/.10/.20 se contrastan además con cuadratura 64/128, fuentes BIE sobremuestreadas por 8 y secciones desplazadas. Estos órdenes no modifican M=24.',
        f"Búsqueda de colocación desde a=0, sin usar a1 como centro: diferencia máxima en a {max(r['placement_difference'] for r in placement):.3e} (.02/.10).", '',
        'La incertidumbre del operador modal y el residual BIE se declaran como errores medidos, no cotas rigurosas del operador continuo. El resultado es evidencia numérica a resolución fija.', '',
        '## Sensibilidad y asintótica', '',
        'La sensibilidad reajusta a y sigma conjuntamente; sus tablas y `internal_physical_certificates.json` conservan también canales complejos firmados, paredes, BIE y decaimiento. Los casos de baja precisión no se convierten en inexistencia.',
        'Los ajustes de a/epsilon y sigma/epsilon² usan modos físicos. Consistencia con el término principal no prueba Big-O ni analiticidad. El ejemplo 5.4 del artículo es internamente discrepante: su sustitución en (5.47) da -beta/2, frente al -beta/12 impreso; se conserva (2.7).', '',
        '## Procedencia', '',
        'Resultados corrected/v3; carpetas embedded_validation y embedded_validation_v2 son legacy/pre-fix y permanecen intactas. Cada ejecución crea una versión nueva; el fingerprint incluye configuración, código científico, Green y harness. Los contrastes rechazan una procedencia incompatible.',
        'Ver run_manifest.json, verification_gate.json, certificates/, publication_claims.csv, independent_validation.json e independent_placement.json.']
    (run/'results_audit.md').write_text('\n'.join(lines)+'\n')
    y.harness.atomic_json(run/'audit_provenance.json',{'scientific_fingerprint':manifest['fingerprint'],
        'verification_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('smoke','full'),default='smoke')
    parser.add_argument('--check-run',type=Path)
    args=parser.parse_args()
    if args.check_run:
        passed=check_results(args.check_run)
        sensitivity=check_internal_sensitivity(args.check_run)
        write_verification_report(args.check_run)
        if not passed or not sensitivity:raise SystemExit('Independent validation did not pass')
        return
    cfg=y.CONFIG
    if args.stage=='smoke':
        cfg=replace(cfg,epsilon_values=(.02,.10,.20),run_internal_convergence_study=False)
    active=y.harness.new_run(y,cfg,purpose=args.stage)
    y.harness.manifest(y,active,'running')
    try:
        y._run_validation(active)
        passed=check_results(active.output_directory)
        sensitivity=check_internal_sensitivity(active.output_directory)
        summaries=list(csv.DictReader((Path(active.output_directory)/'summary.csv').open()))
        all_modes=bool(len(summaries)==len(active.epsilon_values)
                       and all(r['mode_exists_supported']=='True' for r in summaries))
        y.harness.atomic_json(Path(active.output_directory)/'verification_gate.json',
            {'representative_independent_checks_passed':passed,'internal_sensitivity_passed':sensitivity,
             'all_requested_modes_supported':all_modes,'M':24,'stage':args.stage})
        write_verification_report(active.output_directory)
        if not passed or not sensitivity or not all_modes:raise RuntimeError('Independent physics/mathematics gate failed')
    except BaseException as exc:
        y.harness.manifest(y,active,'interrupted' if isinstance(exc,KeyboardInterrupt) else 'failed');raise
    y.harness.manifest(y,active,'completed')

if __name__=='__main__':main()
