"""Read-only cross-check of final published tables; no solver/sweep."""
from pathlib import Path
import csv,json,math
root=Path(__file__).resolve().parents[3]
base=root/'src/beyn_theorems_validations/results';out={}
for name in ('theorem_2_3_iii','theorem_2_3_iv'):
 run=Path(json.loads((base/name/'runs/latest_run.json').read_text())['directory'])
 rows=list(csv.DictReader((run/'summary.csv').open()));cfg=json.loads((run/'run_manifest.json').read_text())['config'];b=cfg['b']
 checks=[]
 for r in rows:
  e=float(r['epsilon']);k=float(r['kb_numerical'])/b;s=float(r['sigma_numerical']);mu=float(r['mu'])
  assert float(r['lambda_1'])<k*k<float(r['lambda_2'])
  assert abs(k*k+s*s-float(r['lambda_2']))<1e-12
  assert abs(float(r['sigma_over_epsilon_squared'])-s/e**2)<1e-10
  assert abs(float(r['asymptotic_coefficient'])-math.pi**3*mu/b**3)<1e-10
  assert r['mode_exists_supported']=='True'
  checks.append({'epsilon':e,'dispersion_error':abs(k*k+s*s-float(r['lambda_2'])),'leading_relative_error':abs(s/e**2-float(r['asymptotic_coefficient']))/float(r['asymptotic_coefficient'])})
 conv=list(csv.DictReader((run/'internal_convergence.csv').open()))
 out[name]={'run':run.name,'cases':checks,'plots':len(list(run.glob('*.png'))),'zero_convergence_changes':sum(float(r['relative_sigma_shift_from_reference'])==0 for r in conv),
  'global_uniqueness_supported':sum(r.get('uniqueness_supported',r.get('full_problem_uniqueness_supported'))=='True' for r in rows)}
run=Path(json.loads((base/'theorem_2_1/runs/latest_run.json').read_text())['directory'])
rows=list(csv.DictReader((run/'paper_a_sweep_M24.csv').open()));critical=list(csv.DictReader((run/'critical_height.csv').open()))
assert len(rows)==20 and all(r['status']=='zero-supported' and r['nonexistence_supported']=='True' and r['diagnostic_complete']=='True' for r in rows)
out['theorem_2_1']={'final_run':run.name,'plots_in_final_run':len(list(run.glob('*.png'))),'absence_supported_cases':len(rows),'critical_brackets':critical,
'baseline_plot_run':'2026-10-07_11-06-56_full_832cf4c8c71d','baseline_plots':97,'baseline_version':'v10; final supplements v15'}
(Path(__file__).parent/'evidence.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
