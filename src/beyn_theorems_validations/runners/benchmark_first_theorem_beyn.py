"""Small warmed assembly benchmark against the frozen pre-optimization lattice."""
from __future__ import annotations
import importlib.util
import json
import statistics
import sys
import time
from dataclasses import replace
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import first_theorem_beyn as f


def uncached_geometry(epsilon, M, b, a, h, harmonic_order):
    t=f.boundary_nodes(M)
    x,y,xp,yp,xpp,ypp=f.circle_geometry(t,epsilon)
    xi=np.stack([x+h,x-h,x,x]); eta=np.stack([y,y,y+h,y-h])
    X=x[None,:,None]-xi[:,None,:]
    direct=y[None,:,None]-eta[:,None,:]
    image=y[None,:,None]+eta[:,None,:]+2*(b+a)
    dx=x[None,:]-xi; dd=y[None,:]-eta; di=y[None,:]+eta+2*(b+a)
    curvature=(xpp*yp-ypp*xp)/(4*f.PI*(xp*xp+yp*yp))
    return X,direct,image,dx,dd,di,xp,yp,curvature,(None,)*4


def main():
    source=Path(__file__).resolve().parents[1]/'tests/reference/lattice_sums_before_optimization.py'
    spec=importlib.util.spec_from_file_location('reference_lattice',source)
    reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(reference)
    current_lattice=f.lattice.lattice_sums;current_geometry=f.assembly_geometry
    cases=[(eps,k,replace(f.CONFIG,harmonic_order=order))
           for eps,k in ((.02,1.5707956179),(.1,.8+.01j),(.2,1.566637+.00001j))
           for order in (20,40)]
    def run(old):
        f.lattice.lattice_sums=reference.lattice_sums if old else current_lattice
        f.assembly_geometry=uncached_geometry if old else current_geometry
        return [f.assemble_matrix(k,eps,24,cfg) for eps,k,cfg in cases]
    try:
        before=run(True);after=run(False)
        errors=[float(np.linalg.norm(a-b)/np.linalg.norm(a)) for a,b in zip(before,after)]
        assert max(errors)==0., errors
        times={}
        for old,label in ((True,'before'),(False,'optimized')):
            samples=[]
            for _ in range(5):
                start=time.perf_counter();run(old);samples.append(time.perf_counter()-start)
            times[label]=statistics.median(samples)
        report={'M':24,'cases':len(cases),'max_relative_matrix_error':max(errors),
                'seconds_per_six_assemblies':times,
                'assembly_speedup':times['before']/times['optimized'],
                'scope':'Warmed real/complex matrices, H20/H40; not an end-to-end sweep estimate'}
    finally:
        f.lattice.lattice_sums=current_lattice;f.assembly_geometry=current_geometry
    config=f.new_run_config(f.CONFIG,"benchmark")
    directory=Path(config.output_directory)
    (directory/'performance_benchmark.json').write_text(json.dumps(report,indent=2)+'\n')
    f.write_run_manifest(directory,config,'benchmark-completed')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
