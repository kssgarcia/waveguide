import pathlib,json,csv,math,hashlib,collections
base=pathlib.Path('src/beyn_theorems_validations');out={}
read=lambda p:list(csv.DictReader(p.open()))
f=lambda r,k:float(r[k]);b=lambda r,k:r[k]=='True'
for t in ['theorem_2_1','theorem_2_3_iii','theorem_2_3_iv']:
 d=pathlib.Path(json.loads((base/'results'/t/'runs/latest_run.json').read_text())['directory']);s=read(d/'summary.csv');p=read(d/'physical_diagnostics.csv');checks=[]
 for r,ph in zip(s,p):
  eps=f(r,'epsilon');kb=f(r,'kb_numerical');sig=f(r,'sigma_numerical');cut=math.pi/2 if t=='theorem_2_1' else math.pi
  assert abs(sig*sig-(cut*cut-kb*kb))<1e-12
  assert abs(sig/eps**2-f(r,'sigma_over_epsilon_squared'))<1e-10
  assert abs(abs(sig-f(r,'sigma_asymptotic'))/abs(f(r,'sigma_asymptotic'))-f(r,'relative_error_sigma'))<1e-12
  assert all(math.isfinite(f(r,k)) for k in ['normalized_residual','equilibrated_residual','equilibrated_sigma_min','final_drop_factor'])
  assert f(r,'normalized_residual')<=1e-8 and f(r,'equilibrated_residual')<=1e-8
  assert f(r,'equilibrated_sigma_min')<=(1e-8 if t=='theorem_2_3_iv' else 1e-6)
  assert f(r,'final_drop_factor')>=100 and f(r,'relative_singular_value_final')<=1e-4
  assert f(ph,'wall_relative_residual')<=1e-6 and f(ph,'boundary_integral_relative_residual')<=5e-3
  assert f(ph,'decay_relative_error_left')<=.2 and f(ph,'decay_relative_error_right')<=.2
  assert b(ph,'monotone_decay_left') and b(ph,'monotone_decay_right')
  if t!='theorem_2_1': assert f(ph,'first_open_channel_relative_amplitude')<=1e-5
  if t=='theorem_2_3_iii':
   assert all(f(ph,k)<=1e-6 for k in ['centerline_relative_residual','odd_parity_relative_residual'])
   assert f(r,'full_vector_residual')<=1e-8 and f(r,'parity_leakage_residual')<=1e-7
   assert r['beyn_estimated_rank']==r['beyn_threshold_rank']==r['beyn_gap_rank']=='1'
   assert not b(r,'full_problem_uniqueness_supported')
  if t=='theorem_2_3_iv':
   assert f(r,'root_residual_norm')<=2e-8 and b(r,'root_success')
   assert all(f(ph,k)<=1e-5 for k in ['boundary_x_even_residual','field_x_even_residual'])
   assert f(ph,'modal_projection_discrepancy')<=1e-7
   assert f(ph,'open_channel_phase_consistency')<=.05 or f(ph,'first_open_channel_relative_amplitude')<=1e-5
   assert abs((f(r,'a_numerical')-eps*f(r,'a1_formula_2_7'))/eps**2-f(r,'scaled_a_correction'))<1e-10
   assert not b(r,'uniqueness_supported')
  else:
   assert f(r,'beyn_final_s0_relative_change')<=1e-3 and f(r,'beyn_final_s1_relative_change')<=1e-3
   assert b(r,'beyn_rank_stable') and b(r,'tighter_margin_rank_consistent')
  assert b(r,'mode_exists_supported')
  checks.append({'epsilon':eps,'spectral_and_physical_recomputed':True})
 metrics={k:max(f(r,k) for r in s) for k in ['normalized_residual','equilibrated_residual','relative_error_sigma']}
 metrics.update({k:max(f(r,k) for r in p) for k in ['wall_relative_residual','boundary_integral_relative_residual','decay_relative_error_left','decay_relative_error_right']})
 if t!='theorem_2_1':metrics['first_open_channel_relative_amplitude']=max(f(r,'first_open_channel_relative_amplitude') for r in p)
 if t!='theorem_2_3_iv':metrics['s0_change']=max(f(r,'beyn_final_s0_relative_change') for r in s);metrics['s1_change']=max(f(r,'beyn_final_s1_relative_change') for r in s)
 info=dict(run=str(d),cases=checks,maxima=metrics,sha256={name:hashlib.sha256((d/name).read_bytes()).hexdigest() for name in ['summary.csv','physical_diagnostics.csv','internal_convergence.csv']})
 if t=='theorem_2_1':
  g=read(d/'paper_a_sweep_M24.csv');info['geometry_counts']=dict(collections.Counter(r['status'] for r in g));absent=[r for r in g if b(r,'nonexistence_supported')]
  assert all(b(r,'diagnostic_complete') and b(r,'beyn_moments_converged') and r['beyn_rank']==r['beyn_threshold_rank']==r['beyn_gap_rank']=='0' and f(r,'contour_solve_residual')<=1e-10 and int(r['local_mode_count'])==0 for r in absent)
  exclusions=[]
  for r in absent:
   a=f(r,'a');ep=f(r,'epsilon');cert=d/'geometry_diagnostics'/f'a_{a:.12g}_epsilon_{ep:.12g}'/'contour_exclusion.json'
   if cert.exists():
    q=json.loads(cert.read_text());assert q['supported'] and q['complete'];assert len(q['history'])>=2 and all(v['resolved'] and v['count']==0 for v in q['history'][-2:]);exclusions.append(str(cert))
  info['absent_cases_recomputed']=len(absent);info['independent_exclusions']=exclusions;info['suspicious_unresolved_absence_cases']=sum(int(r['svd_scan_suspicious_unresolved'])>0 for r in absent)
  info['critical_height']=read(d/'critical_height.csv')
 if t=='theorem_2_3_iv':info['sensitivity_status_counts']=dict(collections.Counter(r['status'] for r in read(d/'internal_convergence.csv')))
 out[t]=info
pathlib.Path('docs/audits/manuscript-reconstruction/recomputed_acceptance.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({t:{k:v for k,v in r.items() if k not in ['cases','sha256','independent_exclusions','critical_height']} for t,r in out.items()},indent=2))
