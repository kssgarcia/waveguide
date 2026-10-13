"""Read-only scientific audit controls for Y; no full sweep."""
from pathlib import Path
import sys,json,csv,math,hashlib
from dataclasses import replace
from unittest.mock import patch
import numpy as np
root=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(root/'src/beyn_theorems_validations'))
import second_theorem_Y as y
if y.ASSEMBLY_VERSION != '2026-10-07-periodic-green-direct-a1-v2':
 raise SystemExit('Este control corresponde a la auditoría legacy/pre-fix. Su evidencia original está en evidence.json; para el código corregido ejecute tests/test_second_theorem_Y_science.py.')
p=Path(__file__).parent
out={}
rows=[]
for beta in (.001,.002,.005,.02,.05,.1,.5):
 row,c,sx,sy=y._solve_reference_mfs(replace(y.CONFIG,shape_beta=beta),160)
 t=2*np.pi*(np.arange(2048)+.173)/2048
 X,Y,Xp,Yp,*_=y.reference_obstacle_geometry(t,replace(y.CONFIG,shape_beta=beta))
 dx=X[:,None]-sx;dy=Y[:,None]-sy;r2=dx*dx+dy*dy;speed=np.hypot(Xp,Yp)
 nx=-Yp/speed;ny=Xp/speed
 off=-(dx*nx[:,None]+dy*ny[:,None])/(2*np.pi*r2)@c
 psi=-np.log(r2)@c/(4*np.pi)
 f=.5*(Y-psi);f1=(f+np.cos(t))/beta
 dt=2*np.pi/len(t);mu_identity=(row.area+dt*np.sum(ny*psi*speed))/(2*np.pi)
 rows.append({'beta':beta,'a1':row.a1_formula_2_7,'a1_over_beta':row.a1_formula_2_7/beta,
              'offgrid_Neumann_residual':float(np.linalg.norm(off-ny)/np.linalg.norm(ny)),
              'mu_identity_relative_error':abs(mu_identity-row.mu)/abs(row.mu),
              'f1_mean':float(np.mean(f1)),'f1_cos2_coefficient':float(2*np.mean(f1*np.cos(2*t)))})
out['small_beta_MFS']=rows
# Direct evaluation of the paper's written f0/f1 in (5.47), to first order.
t=2*np.pi*(np.arange(8192)+.5)/8192;dt=2*np.pi/len(t)
def paper_expansion(beta,constant):
 X=np.sin(t)-beta/2*np.sin(2*t);Y=-np.cos(t)+beta/2*np.cos(2*t)
 Xp=np.cos(t)-beta*np.cos(2*t);Yp=np.sin(t)-beta*np.sin(2*t)
 f=-np.cos(t)+beta*(.5*np.cos(2*t)+constant)
 return dt*np.sum((Y*Xp-3*X*Yp)*f)/(3*np.pi)
beta=1e-5
out['paper_f1_substitution']={'printed_f1_constant':-.25,'a1_slope_from_printed_f1':(paper_expansion(beta,-.25)-paper_expansion(-beta,-.25))/(2*beta),
 'printed_example_slope':-1/12}
import sympy as sp
u,B=sp.symbols('u B',real=True)
Xs=sp.sin(u)-B*sp.sin(2*u)/2;Ys=-sp.cos(u)+B*sp.cos(2*u)/2
fs=-sp.cos(u)+B*(sp.cos(2*u)/2-sp.Rational(1,4))
linear=sp.expand(sp.diff((Ys*sp.diff(Xs,u)-3*Xs*sp.diff(Ys,u))*fs,B).subs(B,0))
out['paper_f1_substitution']['symbolic_first_order_integral']=str(sp.integrate(linear,(u,0,2*sp.pi)))
out['paper_f1_substitution']['symbolic_a1_slope']=str(sp.integrate(linear,(u,0,2*sp.pi))/(3*sp.pi))
# A real failure must not be indistinguishable from an ordinary objective value.
with patch.object(y,'full_singular_metrics',side_effect=ValueError('injected')):
 out['failure_penalties']={'singularity_measure':y.singularity_measure_from_sigma(.2,.1,-.02,y.CONFIG),
                          'relative_measure':y.relative_sv_from_sigma(.2,.1,-.02,y.CONFIG)}
with patch.object(y,'near_zero_eigenpair',side_effect=ValueError('injected')):
 out['failure_penalties']['root_residual']=y._root_residual(np.array([0.,0.]),.1,-.2,.2,y.CONFIG).tolist()
# Every scan value NaN yields no brackets and no extra modes.
with patch.object(y,'build_uniqueness_grid',return_value=np.linspace(1.8,2.8,9)),patch.object(y,'full_singular_metrics',return_value=(math.nan,1.,math.nan,math.nan)):
 mode=type('Mode',(),{'a_numerical':-.02,'kb_numerical':3.})()
 scan,values,candidates=y.whole_band_uniqueness_screen(.1,1.,mode,y.CONFIG)
 out['nan_uniqueness_screen']={'nonfinite':int(np.sum(~np.isfinite(values))),'candidate_count':len(candidates),'validation_pass_formula':sum(c.nonradiating_resolved_candidate for c in candidates)==0}
# Failed least_squares is still selected if its endpoint eigenvalue is tiny.
fake=type('Result',(),{'success':False,'x':np.array([0.,0.]),'nfev':80})()
with patch.object(y,'least_squares',return_value=fake),patch.object(y,'near_zero_eigenpair',return_value=(1e-12+1e-12j,np.ones(24),1.)):
 result,source,lam,norm=y._solve_complex_eigenvalue_root(.1,-.2,.2,2.,[('theorem-leading',np.array([0.,0.]))],y.CONFIG)
 out['failed_root_control']={'optimizer_success':result.success,'selected':source,'root_norm':norm,'root_success_flag_formula':norm<=y.CONFIG.root_residual_tolerance}
base=root/'src/beyn_theorems_validations/results/theorem_2_3_iv/embedded_validation'
saved=list(csv.DictReader((base/'summary.csv').open()));physical=list(csv.DictReader((base/'physical_diagnostics.csv').open()))
out['saved_results']=[{'epsilon':float(r['epsilon']),'a':float(r['a_numerical']),'sigma':float(r['sigma_numerical']),'expected_bic_resolved':r['expected_bic_resolved'],'physical_bic_verified':r['physical_bic_verified'],'scan_run':r['uniqueness_screen_run'],'scan_passed':r['uniqueness_screen_passed']} for r in saved]
out['saved_physical']=[{'epsilon':float(r['epsilon']),'open_channel_ratio':float(r['first_open_channel_relative_amplitude']),
                       'BIE_residual':float(r['boundary_integral_relative_residual']),'decay_relative_error':max(float(r['decay_relative_error_left']),float(r['decay_relative_error_right']))} for r in physical]
paths=['src/beyn_theorems_validations/second_theorem_Y.py','src/beyn_theorems_validations/execution_harness.py',
       'docs/references/WG_03.03.25.pdf','docs/references/linton1998greens.pdf',
       'src/beyn_theorems_validations/results/theorem_2_3_iv/embedded_validation/summary.csv',
       'src/beyn_theorems_validations/results/theorem_2_3_iv/embedded_validation/physical_diagnostics.csv']
out['sha256']={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in paths}
(p/'evidence.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
