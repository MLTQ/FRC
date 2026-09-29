"""Frozen axisymmetric manufactured moving-current rate protocol."""
from dataclasses import asdict,replace
from pathlib import Path
from time import perf_counter
import argparse
import hashlib
import numpy as np
from scipy.constants import electron_mass as ME
from plasma_bench.moving_current import MovingCurrent
from .moving_species_reference import MovingSpeciesReference
from .moving_current_case import initialize,norms
from .moving_current_records import provenance
from .moving_current_work import verify as verify_work
from .verification_records import write_record


def protocol():
    ref=MovingSpeciesReference();cases=[]
    def add(name,nr,reference=ref,**changes):
        case=dict(name=name,nr=nr,reference=asdict(reference),phase=.23,quadrature=5,tolerance=1e-10,max_iterations=2000)
        case.update(changes);cases.append(case)
    for nr in (8,16,32,64,128):add(f'moving_{nr}',nr)
    for nr in (8,16,32,64):add(f'target_{nr}',nr,replace(ref,density=1e20))
    for nr in (8,16,32,64):add(f'massless_{nr}',nr,replace(ref,electron_mass=0.))
    add('pair_16',16,replace(ref,ion_mass=ME))
    for phase in (0.,1.17):add(f'phase_{phase:g}',16,phase=phase)
    add('quadrature_16',16,quadrature=7);add('tolerance_16',16,tolerance=1e-12)
    add('tolerance_64',64,tolerance=2e-12)
    return cases


def measure(case):
    ref,model,psi,g,ue,ui,v,nt,source,time=initialize(case);geo=model.geometry
    solver=MovingCurrent(model,tolerance=case['tolerance'],max_iterations=case['max_iterations'])
    arrays=dict(psi=psi,g=g,number_density=model.electron.density,electron_entropy=ue,ion_entropy=ui,
        regular_velocity=v,density_rate=nt,radial_faces=model.mesh.rf,axial_faces=model.mesh.zf)
    original={k:a.copy() for k,a in arrays.items()}
    try:rate=solver.rate(psi,g,ue,ui,v,nt,electric_source=source)
    except RuntimeError as error:
        assert all(np.array_equal(a,arrays[k]) for k,a in original.items())
        return dict(completed=False,status='rejected',failure=str(error),input_unchanged=True),arrays
    p,q=geo.unpack(rate.induction.magnetic)
    arrays.update(rate_psi=p,rate_g=q,rate_current=rate.induction.current,electric_full=rate.electric_full,
        electric_reduced=rate.induction.electric,potential=rate.potential)
    weak,_=solver.terms(psi,g,ue,ui,v,nt);forcing=sum(weak.values())+model.mass@geo.project_electric(source)
    work=verify_work(model,psi,g,ue,ui,nt,rate.induction,forcing)
    metrics=norms(model.mesh,ref,time,p,q,rate.induction.current,rate.electric_full)
    omitted={}
    if case['name']=='moving_32':
        for name in ('density','vorticity','bulk'):
            result=solver.solver.solve(forcing-weak[name]);change=result.magnetic-rate.induction.magnetic
            difference=np.sqrt(change@(model.hessian@change)/(rate.induction.magnetic@(model.hessian@rate.induction.magnetic)))
            omitted[name]=dict(magnetic_rate_relative_change=float(difference),solver_residual=result.residual)
            arrays['without_'+name+'_rate']=result.magnetic
    assert all(np.array_equal(a,arrays[k]) for k,a in original.items())
    return dict(completed=True,status='succeeded',failure=None,input_unchanged=True,
        iterations=rate.induction.iterations,linear_residual=rate.induction.residual,
        linear_refinements=rate.induction.refinements,residual_components=rate.induction.residual_components,metrics=metrics,
        work=work,omission_controls=omitted),arrays


def run(directory,group=None):
    directory.mkdir(parents=True,exist_ok=True);source=provenance()
    if group is not None and group not in {c['name'] for c in protocol()}:raise ValueError('Unknown case')
    for case in protocol():
        if group is not None and case['name']!=group:continue
        start=perf_counter();record,arrays=measure(case);assert source==provenance()
        path=directory/(case['name']+'.npz');np.savez_compressed(path,**arrays)
        record.update(configuration=case,provenance=source,elapsed_seconds=perf_counter()-start,
            fields=dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        write_record(directory/(case['name']+'.json'),record)
        print(case['name'],record['status'],record.get('metrics'),record['failure'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',default='outputs/moving_current_verification');parser.add_argument('--group')
    args=parser.parse_args();run(Path(args.directory),args.group)
