"""Scientific gate then the full X sweep, sharing signed case checkpoints."""
import sys,json,shutil
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import second_theorem_X as x

def run(config):
    output=Path(config.output_directory)
    mu,nu,rows,circle=x.run_mu_convergence(config)
    checks=[]
    for epsilon in (.02,.10,.20):
        result,*_=x.harness.case(x,config,f'epsilon_{epsilon:.8f}',(epsilon,mu,nu))
        passed=bool(result.mode_exists_supported and result.odd_uniqueness_supported
                    and result.diagnostic_complete and result.physical_embedded_mode_verified)
        checks.append({'epsilon':epsilon,'passed':passed,'result':x.asdict(result)})
        x.harness.atomic_json(output/'representative_gate.json',checks)
        print(f'gate epsilon={epsilon}: {passed}',flush=True)
        if not passed:raise RuntimeError(f'X scientific gate failed at epsilon={epsilon}; see certificate')
    x._run_validation(config)
    summary=output/'summary.csv'
    import csv
    records=list(csv.DictReader(summary.open()))
    for r in records:
        if r['mode_exists_supported']!='True' or r['odd_uniqueness_supported']!='True':
            raise RuntimeError(f'Incomplete corrected result at epsilon={r["epsilon"]}')
    x.harness.atomic_json(output/'scientific_audit.json',{
        'assembly_version':x.ASSEMBLY_VERSION,'M':config.bem_M,'cases':len(records),
        'existence_and_odd_uniqueness_passed':True,'full_problem_uniqueness_established':False,
        'scope':'Physical BIC and numerical odd-sector uniqueness inside sampled contours; even sector and omitted thresholds not excluded.'})

if __name__=='__main__':x.harness.run_entrypoint(x,x.CONFIG,run)
