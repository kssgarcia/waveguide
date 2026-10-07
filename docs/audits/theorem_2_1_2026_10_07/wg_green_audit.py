from pathlib import Path
import sys, json, time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[3] / 'src/beyn_theorems_validations'))
import lattice_sums as l
start=time.time(); out=[]
# Independent Fourier representation of the periodic Green function (Linton 2.13).
def periodic(X,Y,k,d,beta=0.,N=800):
 m=np.arange(-N,N+1); bm=beta+2*np.pi*m/d
 g=np.sqrt(bm**2-k*k+0j); opened=np.abs(bm)<k.real
 g[opened]=-1j*np.sqrt(k*k-bm[opened]**2+0j)
 return -np.sum(np.exp(-g*abs(X)+1j*bm*Y)/g)/(2*d)
def modal(X,y,eta,k,d=2.,N=800):
 m=np.arange(1,N+1); q=np.sqrt((m*np.pi/d)**2-k*k+0j)
 return -np.sum(np.exp(-q*abs(X))*np.sin(m*np.pi*y/d)*np.sin(m*np.pi*eta/d)/q)/d
for k in [1.,1.+.02j,1.-.02j,1.5704,1.5704+.0001j]:
 for H in [10,20,40]:
  S=l.lattice_sums(4.,k,0.,200,H)
  for X,y,eta in [(.2,1.65,1.55),(.4,.7,1.6),(1.5,.2,1.7)]:
   v=l.greens_dirichlet(X,y,0.,eta,S,k,2.)
   ref=modal(X,y,eta,complex(k))
   out.append(dict(k=str(k),H=H,X=X,y=y,eta=eta,error_abs=abs(v-ref),error_rel=abs(v-ref)/abs(ref)))
S=l.lattice_sums(1.,2.,np.sqrt(2),200,20)
v=l.greens_periodic(0.,.01,S,2.,1.)
print('Linton table2',v,'error',abs(v-(-.4595298795-.3509130869j)))
for H in [20,40]:
 rows=[r for r in out if r['H']==H]
 print('H',H,'max abs',max(r['error_abs'] for r in rows),'max rel',max(r['error_rel'] for r in rows))
for r in out:
 if r['error_rel']>1e-4: print(r)
json.dump(out,open('/tmp/wg_green_audit.json','w'),indent=2)
print('elapsed',time.time()-start)
