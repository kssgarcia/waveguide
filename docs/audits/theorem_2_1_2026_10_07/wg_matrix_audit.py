from pathlib import Path
exec(open(Path(__file__).with_name('wg_green_audit.py')).read().split('for k in [')[0])
import first_theorem_beyn as f
from dataclasses import replace
from scipy.special import hankel1
import csv
# Accelerated independent Dirichlet modal series: subtract the k=0 tail analytically.
def gm(x,y,xi,eta,k,N=3000):
 d=2.; X=abs(x-xi); m=np.arange(1,N+1); p=m*np.pi/d; q=np.sqrt(p*p-k*k+0j)
 delta=(np.exp(-q*X)/q-np.exp(-p*X)/p)/d
 smooth=-np.sum(delta*np.sin(p*y)*np.sin(p*eta))
 r=np.exp(-np.pi*X/d)
 def dist(z): return (1-r)**2+4*r*np.sin(np.pi*z/(2*d))**2
 return smooth-np.log(dist(y+eta)/dist(y-eta))/(4*np.pi)
def independent(e,k,M=24,N=3000):
 t=f.boundary_nodes(M); A=np.eye(M,dtype=complex); h=1e-6
 for i,p in enumerate(t):
  x,y=e*np.cos(p),e*np.sin(p)+1.6
  for j,s in enumerate(t):
   xi,eta=e*np.cos(s),e*np.sin(s)+1.6
   def g(xx,yy):
    v=gm(x,y,xx,yy,k,N)
    if i==j: v-= -1j/4*hankel1(0,k*np.hypot(x-xx,y-yy))
    return v
   gx=(g(xi+h,eta)-g(xi-h,eta))/(2*h)
   gy=(g(xi,eta+h)-g(xi,eta-h))/(2*h)
   K=-e*np.sin(s)*gy-e*np.cos(s)*gx
   if i==j: K+=-1/(4*np.pi)
   A[i,j]-=4*np.pi/M*K
 return A
rows=list(csv.DictReader(open(str(Path(__file__).resolve().parents[3] / 'src/beyn_theorems_validations/theorem_2_1_publication_ready/summary.csv'))))
for e in [.02,.1,.2]:
 r=next(x for x in rows if float(x['epsilon'])==e); k=float(r['kb_numerical'])
 B=independent(e,k); svB=np.linalg.svd(B,compute_uv=False)
 for H in [20,40]:
  A=f.assemble_matrix(k,e,24,replace(f.CONFIG,harmonic_order=H))
  u,s,v=np.linalg.svd(A); vec=v.conj().T[:,-1]
  print('e',e,'H',H,'matrix_relative',np.linalg.norm(A-B)/np.linalg.norm(B),'max_entry',np.max(abs(A-B)),'sv_min',s[-1],'modal_svmin',svB[-1],'modal_vector_residual',np.linalg.norm(B@vec)/np.linalg.norm(B))
print('elapsed',time.time()-start)
