#!/usr/bin/env python3
"""Reconstruct one certified discrete mode at epsilon=.10,a=.60; no root search."""
from pathlib import Path
from datetime import datetime,timezone
import sys,json
import numpy as np
import matplotlib.pyplot as plt
import generate_publication as pub

pub.check_sources()
sys.path.insert(0,str(pub.ROOT/'src/beyn_theorems_validations'))
import first_theorem_beyn as solver
import second_theorem_Y as batched
cfg=pub.read_config(solver,'theorem_2_1')
r=pub.csv('theorem_2_1','summary.csv');r=r[np.isclose(r.epsilon,.1)].iloc[0]
eps=float(r.epsilon);a=float(r.a);kb=float(r.kb_numerical);sigma=float(r.sigma_numerical)
assert bool(r.mode_exists_supported)
t=solver.boundary_nodes(cfg.bem_M)
A=solver.assemble_matrix(kb,eps,cfg.bem_M,cfg)
_,sv,vh=np.linalg.svd(A);v=vh[-1].conj();v/=np.linalg.norm(v);v*=np.exp(-1j*np.angle(v[np.argmax(abs(v))]))
res=float(np.linalg.norm(A@v)/(np.linalg.norm(A,2)*np.linalg.norm(v)))
assert res<1e-8
x=np.linspace(-4,4,121);y=np.linspace(-cfg.b,cfg.b,71);xx,yy=np.meshgrid(x,y);pts=np.column_stack([xx.ravel(),yy.ravel()])
geom=(eps*np.cos(t),a+eps*np.sin(t),-eps*np.sin(t),eps*np.cos(t))
mask=np.hypot(pts[:,0],pts[:,1]-a)<eps+.025*cfg.b
field=np.full(len(pts),np.nan+1j*np.nan,complex);valid=np.flatnonzero(~mask)
far=valid[np.min(abs(pts[valid,0,None]-geom[0]),axis=1)>=cfg.b];near=np.setdiff1d(valid,far)
for chunk in np.array_split(far,max(1,int(np.ceil(len(far)/180)))):
 if len(chunk):field[chunk]=pub.modal(pts[chunk],v,geom,kb,cfg.b)
for chunk in np.array_split(near,max(1,int(np.ceil(len(near)/80)))):
 if len(chunk):field[chunk]=2*np.pi/len(t)*(batched._offgrid_kernel_array(pts[chunk,0],pts[chunk,1],*geom,kb,cfg)@v)
# Original discrete solver uses y relative to the displaced circle.
picks=valid[np.linspace(0,len(valid)-1,7,dtype=int)];G=solver.make_extended_field_green(kb,cfg)
direct=solver.reconstruct_field_at_points([(float(pts[i,0]),float(pts[i,1]-a)) for i in picks],v,t,eps,cfg,G)
discrepancy=float(np.linalg.norm(direct-field[picks])/max(np.linalg.norm(direct),1e-30));assert discrepancy<1e-7
yn,yw=np.polynomial.legendre.leggauss(96);yn*=cfg.b;yw*=cfg.b
dist=np.unique(np.r_[1.5,2,2.5,3,1.5+np.linspace(0,4.5,24)/sigma]);sections=np.array([[pub.modal(np.column_stack([np.full(len(yn),sign*d),yn]),v,geom,kb,cfg.b) for d in dist] for sign in [-1,1]])
norms=np.sqrt(np.sum(abs(sections)**2*yw,axis=2));assert np.all(np.diff(norms,axis=1)<0)
ix=np.flatnonzero(np.isin(dist,[1.5,2,2.5,3]));rates=[-np.polyfit(dist[ix],np.log(n[ix]),1)[0] for n in norms];assert max(abs(np.array(rates)-sigma))/sigma<1e-3
folder=pub.PAPER/'data/field_reconstructions'/('theorem_2_1_'+datetime.now(timezone.utc).strftime('%Y-%m-%d_%H-%M-%S'));folder.mkdir(parents=True,exist_ok=False)
np.savez_compressed(folder/'field.npz',x=x,y=y,field=field.reshape(xx.shape),mask=mask.reshape(xx.shape),theta=t,vector=v,distances=dist,section_y=yn,section_weights=yw,section_fields=sections,section_norms=norms)
meta={'epsilon':eps,'a':a,'kb':kb,'sigma':sigma,'M':cfg.bem_M,'normalized_residual':res,'scalar_batched_discrepancy':discrepancy,'measured_decay_rates':rates,'near_boundary_mask_width_over_b':.025,'root_source':str(pub.RUNS['theorem_2_1']/'summary.csv'),'source_lock_sha256':pub.sha(pub.PAPER/'scripts/source_lock.json'),'data_sha256':pub.sha(folder/'field.npz'),'method':'Full matrix SVD at saved root; centered-difference Green near field and analytic modal far field; original scalar implementation independently sampled; no root refinement.'}
fig,axs=plt.subplots(2,1,figsize=(7.1,5.2),layout='constrained');z=field.reshape(xx.shape);scale=np.nanmax(abs(z));im=axs[0].pcolormesh(x,y,np.ma.array(z.real/scale,mask=mask.reshape(xx.shape)),cmap='RdBu_r',vmin=-1,vmax=1,shading='nearest',rasterized=True)
axs[0].add_patch(plt.Circle((0,a),eps,color='#777777'));axs[0].set(xlabel='$x/b$',ylabel='$y/b$',title=fr'Theorem 2.1: $\varepsilon={eps:g}$, $a/b={a:g}$');axs[0].set_aspect('equal');fig.colorbar(im,ax=axs[0],label=r'$\Re(u)/\max|u|$',shrink=.75)
for i in range(2):axs[1].semilogy(sigma*(dist-dist[0]),norms[i]/norms[i,0],marker=['o','s'][i],ms=3,ls=['-','--'][i],label=['Left sections','Right sections'][i])
axs[1].semilogy(sigma*(dist-dist[0]),np.exp(-sigma*(dist-dist[0])),color='black',ls=':',label=r'$e^{-\sigma(|x|-x_0)}$');axs[1].set(xlabel=r'$\sigma(|x|-x_0)$',ylabel='Normalized transverse L2 norm');axs[1].legend()
for ext in ['pdf','png']:fig.savefig(pub.PAPER/'figures'/f'field_theorem_2_1.{ext}',dpi=300)
plt.close(fig);meta['figure_sha256']={ext:pub.sha(pub.PAPER/'figures'/f'field_theorem_2_1.{ext}') for ext in ['pdf','png']};(folder/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
