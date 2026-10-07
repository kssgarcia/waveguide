from pathlib import Path
import sys,io,contextlib
import numpy as np
from dataclasses import replace
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[3] / 'src/beyn_theorems_validations'))
import first_theorem_beyn as f
cfg=replace(f.CONFIG,bem_M=4,beyn_probe_dim=3)
for roots in [[],[.5],[.5,.8]]:
 def mat(z,e,M,c):
  a=np.ones(M,dtype=complex)
  for j,r in enumerate(roots): a[j]=z-r
  return np.diag(a)
 with patch.object(f,'assemble_matrix',mat),contextlib.redirect_stdout(io.StringIO()):
  raw,_,diag,_,_=f.beyn_discover(.1,192,cfg)
 print('toy roots',roots,'rank',diag.estimated_rank,'recovered',raw)
for case in ['all_nan','failed_refinement']:
 def objective(s,*args):
  return np.nan if case=='all_nan' else 1+((np.log(s)-np.log(.02))**2)
 with patch.object(f,'absolute_singular_value_from_sigma',objective),patch.object(f,'run_candidate_refinement',side_effect=RuntimeError('injected failure')),contextlib.redirect_stdout(io.StringIO()):
  scan,_=f.independent_whole_band_svd_scan(.1,.35,f.CONFIG,'fault-injection')
 print(case,'minima',scan.local_minima_count,'refined',scan.refined_candidate_count,'nonexistence_supported',scan.nonexistence_supported,'status',scan.status)
# Rank gap can hide a genuine second residue, unlike the empty-contour condition.
print('rank test',f.estimate_beyn_rank(np.array([1.,1e-4,1e-10]),f.CONFIG))
print('skewed residue rank',f.estimate_beyn_rank(np.array([1.,1e-5,1e-7,1e-9]),f.CONFIG))
def ab(s,*args): return 1+1e8*(s-.02)**2
def metrics(k,e,M,c):
 s=f.sigma_from_kb(k,c);v=ab(s);return v,1e10,v/1e10
with patch.object(f,'absolute_singular_value_from_sigma',ab),patch.object(f,'singular_metrics',metrics),contextlib.redirect_stdout(io.StringIO()):
 _,r=f.run_candidate_refinement(.1,1,f.DiscoverySeed(f.kb_from_sigma(.02,f.CONFIG),0,'mock-nonsingular'),(.01,.03),f.CONFIG)
print('everywhere nonsingular accepted',r.resolved,'svmin',r.sigma_min_final,'relative',r.relative_singular_value_final,'drop',r.final_drop_factor)
