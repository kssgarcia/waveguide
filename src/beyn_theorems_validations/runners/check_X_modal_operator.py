"""Independent accelerated modal operator check at corrected roots; no sweep."""
import sys,argparse,json
from pathlib import Path
from dataclasses import replace
import numpy as np
from scipy.special import hankel1
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import second_theorem_X as x


def modal_matrix(kb,epsilon,config,terms=3000):
    """Linton modal Green, k=0 tail summed analytically; separate FD kernel."""
    width=2*config.b;k=complex(kb)/config.b
    n=np.arange(1,terms+1);p=n*np.pi/width
    q=np.where(p<k.real,-1j*np.sqrt(k*k-p*p+0j),np.sqrt(p*p-k*k+0j))
    def green(xx,y,xi,eta,regularized=False):
        X=abs(xx-xi);yf=y+config.b;ys=eta+config.b
        correction=(np.exp(-q*X)/q-np.exp(-p*X)/p)/width
        smooth=-np.sum(correction*np.sin(p*yf)*np.sin(p*ys))
        r=np.exp(-np.pi*X/width)
        one_minus_r=-np.expm1(-np.pi*X/width)
        def distance(z):return one_minus_r**2+4*r*np.sin(np.pi*z/(2*width))**2
        value=smooth-np.log(distance(y+eta+2*config.b)/distance(y-eta))/(4*np.pi)
        if regularized:value+=.25j*hankel1(0,k*np.hypot(xx-xi,y-eta))
        return value
    t=x.boundary_nodes(config.bem_M);xx,yy,xp,yp,xpp,ypp=x.obstacle_geometry(t,epsilon,config)
    A=np.eye(config.bem_M,dtype=complex);h=config.finite_difference_step
    for i in range(len(t)):
        for j in range(len(t)):
            reg=i==j
            gx=(green(xx[i],yy[i],xx[j]+h,yy[j],reg)-green(xx[i],yy[i],xx[j]-h,yy[j],reg))/(2*h)
            gy=(green(xx[i],yy[i],xx[j],yy[j]+h,reg)-green(xx[i],yy[i],xx[j],yy[j]-h,reg))/(2*h)
            weighted=xp[j]*gy-yp[j]*gx
            if reg:weighted+=(xpp[j]*yp[j]-ypp[j]*xp[j])/(4*np.pi*(xp[j]**2+yp[j]**2))
            A[i,j]-=4*np.pi/len(t)*weighted
    Q=x.odd_lift_matrix(config.bem_M)
    return Q.conj().T@A@Q


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('run',type=Path)
    args=parser.parse_args();run=args.run.resolve()
    manifest=json.loads((run/'run_manifest.json').read_text())
    cfg=x.Config(**{k:tuple(v) if isinstance(v,list) else v for k,v in manifest['config'].items()})
    results=[]
    for epsilon in (.02,.1,.2):
        cert=json.loads((run/'geometry_diagnostics'/f'epsilon_{epsilon:.12g}'/'certificate.json').read_text())
        kb=cert['kb_numerical'];_,_,_,v,_=x.odd_singular_pair(kb,epsilon,24,cfg)
        A=x.assemble_odd_matrix(kb,epsilon,24,cfg)
        B=modal_matrix(kb,epsilon,cfg,12000);C=modal_matrix(kb,epsilon,cfg,24000)
        norm=np.linalg.norm(C,2)
        result={'epsilon':epsilon,'kb':kb,'M':24,'modal_terms':[12000,24000],
                'matrix_relative_difference':float(np.linalg.norm(A-C)/np.linalg.norm(C)),
                'modal_truncation_relative_change':float(np.linalg.norm(B-C)/np.linalg.norm(C)),
                'modal_vector_residual':float(np.linalg.norm(C@v)/(norm*np.linalg.norm(v)))}
        results.append(result)
        x.harness.atomic_json(run/'independent_modal_operator.json',results)
        print(result,flush=True)

if __name__=='__main__':main()
