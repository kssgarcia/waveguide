import sys,json,csv,math
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch
import numpy as np
root=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(root/'src/beyn_theorems_validations'))
import second_theorem_X as x
out={}
c=x.CONFIG
with patch.object(x,'refine_candidate_for_M',return_value=(2.,.02,1.,1e10,1e-10,1e4,True)):
 _,mode=x.run_candidate_refinement(.1,1,x.DiscoverySeed(2.,0.,'control'),(.01,.03),c)
 out['nonsingular_valley_control']={'smin':mode.sigma_min_final,'smax':mode.sigma_max_final,'resolved':mode.resolved}
seed=x.DiscoverySeed(2.,0.,'control')
with patch.object(x,'absolute_sv_from_sigma',side_effect=lambda s,*args:1+(s-.02)**2),patch.object(x,'odd_singular_metrics',return_value=(1.,1e10,1e-10)),patch.object(x,'minimize_scalar',return_value=type('R',(),{'success':False,'x':.02,'message':'injected failure'})()):
 _,mode=x.run_candidate_refinement(.1,1,seed,(.01,.03),replace(c,minimum_drop_factor=1.))
 out['failed_optimizer_control']={'success':False,'resolved':mode.resolved}

def modal(X,y,eta,k,N=1200):
 n=np.arange(1,N+1);p=n*np.pi/2;q=np.sqrt(p*p-complex(k)**2+0j)
 q=np.where((p<complex(k).real),-1j*np.sqrt(complex(k)**2-p*p+0j),q)
 return -np.sum(np.exp(-q*abs(X))*np.sin(p*y)*np.sin(p*eta)/q)/2
comparisons=[]
for kb in (2.,2.8,complex(2.8,.01),complex(2.8,-.01)):
 for order in (20,40):
  cfg=replace(c,harmonic_order=order)
  g,_=x.make_green_functions(kb,cfg);gf=x.make_extended_field_green(kb,cfg)
  for xx,y,xi,eta in ((.1,.16,-.1,.19),(.1,.16,-.1,-.19),(.4,.8,-.1,.75)):
   ref=modal(xx-xi,y+1,eta+1,kb)
   comparisons.append({'kb':str(kb),'H':order,'points':[xx,y,xi,eta],'assembly_relative_error':abs(g(xx,y,xi,eta)-ref)/max(abs(ref),1e-30),'field_relative_error':abs(gf(xx,y,xi,eta)-ref)/max(abs(ref),1e-30),'assembly_vs_field':abs(g(xx,y,xi,eta)-gf(xx,y,xi,eta))})
out['green_modal_comparisons']=comparisons
A=x.assemble_full_matrix(2.8,.2,24,c);Q=x.odd_lift_matrix(24);R=x.reflection_matrix(24)
B=Q.conj().T@A@Q
out['odd_projection']={'orthonormal_error':float(np.linalg.norm(Q.conj().T@Q-np.eye(12))),'parity_basis_error':float(np.linalg.norm(R@Q+Q)),'commutator':float(np.linalg.norm(A@R-R@A)/np.linalg.norm(A)),'imaginary_reduced_relative':float(np.linalg.norm(B.imag)/np.linalg.norm(B))}
G=x.make_extended_field_green(2.8,c)
try:out['far_field_limit']=str(G(3.,-.95,-.25,.19))
except Exception as e:out['far_field_limit']={'error':f'{type(e).__name__}: {e}'}
p=root/'src/beyn_theorems_validations/results/theorem_2_3_iii/embedded_validation'
rows=list(csv.DictReader((p/'summary.csv').open()));diags=list(csv.DictReader((p/'beyn_diagnostics.csv').open()))
checks=[]
for r in rows:
 e=float(r['epsilon']);N=int(r['beyn_final_quadrature_points']);d=next(d for d in diags if float(d['epsilon'])==e and int(d['quadrature_points'])==N)
 checks.append({'epsilon':e,'rank':int(d['estimated_rank']),'threshold_rank':int(d['threshold_rank']),'gap_rank':int(d['gap_rank']),'s0_change':float(r['beyn_final_s0_relative_change']),'unique_claim':r['unique_embedded_mode_supported'],'physical':r['physical_embedded_mode_verified'],'relative_sigma_error':float(r['relative_error_sigma']),'equation_error':abs(float(r['k2_numerical'])+float(r['sigma_numerical'])**2-float(r['lambda_2']))})
out['saved_results_checks']=checks
(root/'docs/audits/theorem_2_3_iii_2026_10_13/evidence.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
