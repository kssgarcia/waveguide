from pathlib import Path
import sys,io,contextlib,time,csv,json
import numpy as np
from unittest.mock import patch
from dataclasses import replace
sys.path.insert(0,str(Path(__file__).resolve().parents[3] / 'src/beyn_theorems_validations'))
import first_theorem_beyn as f
l=f.lattice
original=f.make_green_functions
def wrapped(kb,c):
 d=4*c.b; k=kb/c.b; S=l.lattice_sums(d,k,0.,c.lattice_terms,c.harmonic_order)
 def wrap(z): return (z+d/2)%d-d/2
 def g(x,y,xi,eta):
  return l.greens_periodic(x-xi,wrap(y-eta),S,k,d)-l.greens_periodic(x-xi,wrap(y+eta+2*(c.b+c.a)),S,k,d)
 def reg(x,y,xi,eta):
  return l.greens_periodic_reg(x-xi,wrap(y-eta),S,k,d)-l.greens_periodic(x-xi,wrap(y+eta+2*(c.b+c.a)),S,k,d)
 return g,reg
rows=list(csv.DictReader(open(str(Path(__file__).resolve().parents[3] / 'src/beyn_theorems_validations/theorem_2_1_publication_ready/summary.csv')))); out=[]; start=time.time()
for e in [.02,.1,.2]:
 s=float(next(r for r in rows if float(r['epsilon'])==e)['sigma_numerical'])
 for mode in ['original','wrapped']:
  for H in ([20] if mode=='original' else [20,40]):
   c=replace(f.CONFIG,harmonic_order=H)
   with patch.object(f,'make_green_functions',original if mode=='original' else wrapped):
    v=f.refine_candidate_for_M(e,24,.8*s,1.2*s,c)
   r=dict(e=e,mode=mode,H=H,sigma=v[1],svmin=v[2],relative_sv=v[4],interior=v[-1]);out.append(r);print(r,flush=True)
json.dump(out,open('/tmp/wg_wrap_root_audit.json','w'),indent=2)
print('elapsed',time.time()-start)
