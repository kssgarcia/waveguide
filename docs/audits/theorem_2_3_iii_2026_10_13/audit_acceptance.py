from pathlib import Path
exec(Path(__file__).with_name('audit_checks.py').read_text().split('out={}')[0])
p=root/'docs/audits/theorem_2_3_iii_2026_10_13/evidence.json';out=json.loads(p.read_text())
with patch.object(x,'absolute_sv_from_sigma',side_effect=lambda s,*args:1+(s-.02)**2),patch.object(x,'odd_singular_metrics',return_value=(1e-12,1.,1e-12)),patch.object(x,'minimize_scalar',return_value=type('R',(),{'success':False,'x':.02,'message':'injected failure'})()):
 _,mode=x.run_candidate_refinement(.1,1,x.DiscoverySeed(2.,0.,'control'),(.01,.03),x.CONFIG)
 out['failed_optimizer_control']={'success':False,'resolved':mode.resolved,'configuration':'original default tolerances'}
cfg=x.CONFIG;diag=x.BeynDiagnostics(.1,24,12,384,1,2,1,1e6,1.,1e-14,1,1,0)
conv=[x.BeynConvergenceRow(.1,N,1,1.,1e-10,1e-12,.8,1,2.,0.,.1,True) for N in (192,384)]
mode=x.ModeResult(.1,1,'mock',2.,2.,1e-12,1.,1e-12,1e4,math.nan,True)
physical=x.PhysicalDiagnostics(**{field.name:(False if field.type=='bool' else .1) for field in __import__('dataclasses').fields(x.PhysicalDiagnostics)})
with patch.object(x,'run_beyn_convergence_study',return_value=(np.array([2.]),np.array([1.]),diag,[],[x.DiscoverySeed(2.,0.,'mock')],[diag],[],conv,math.nan)),patch.object(x,'run_candidate_refinement',return_value=([],mode)),patch.object(x,'beyn_discover',return_value=(np.array([2.]),np.array([1.]),diag,[],np.eye(12))),patch.object(x,'parity_commutator_residual',return_value=.3),patch.object(x,'physical_mode_diagnostics',return_value=physical):
 result,*_=x.validate_epsilon(.1,1.,0.,cfg)
 out['failed_physical_and_unconverged_count_control']={'unique_embedded_mode_supported':result.unique_embedded_mode_supported,'physical_embedded_mode_verified':result.physical_embedded_mode_verified,'s0_change':result.beyn_final_s0_relative_change,'threshold_rank':2,'gap_rank':1,'parity_commutator':.3}
rows=[]
for rho in (.5,.65,.8):
 for N in (80,160):
  rows.append(__import__('dataclasses').asdict(x.compute_mu_mfs(replace(cfg,mu_mfs_source_scale=rho),N)))
out['mu_source_scale_control']=rows
p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('failed_optimizer_control','failed_physical_and_unconverged_count_control','mu_source_scale_control')},indent=2))
