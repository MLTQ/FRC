"""Mixed radial/axial manufactured and unforced evolution; no plasma adoption."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.cylindrical_entropy_fluid import CylindricalEntropyFluid
from .mixed_euler_reference import MixedEulerReference
from .mixed_euler_step import step
from .mixed_euler_sampling import sample,grid,differences
from .cylindrical_entropy_reference import physical
from .mixed_euler_records import provenance,freeze,save


def initialize(spec):
    ref=MixedEulerReference();d=CylindricalEntropyFluid(ref.radius,2*ref.half_length,spec['elements'],spec['points'],ref.mass,degree=spec.get('degree',5))
    return ref,d,ref.state(d.r,d.z)


def totals(driver,state,ref,t,work,sources):
    f=physical(state,driver.r,ref)
    return np.r_[t,np.sum(driver.volume*f['energy']),np.sum(driver.volume[...,None]*state,axis=(0,1,2)),work,sources]


def protocol():
    specs=[dict(manufactured=mode,elements=n,points=8*n) for mode in (True,False) for n in (2,4,8,16)]
    specs += [dict(manufactured=mode,elements=4,points=32,step_factor=s) for mode in (True,False) for s in (2,4)]
    specs += [dict(manufactured=True,elements=8,points=32),dict(manufactured=True,elements=4,points=64)]
    return specs


def evolve(directory,index,spec,identity):
    directory=Path(directory)
    if (directory/f'fluid_{index:02d}.json').exists():raise ValueError('Terminal result already exists')
    start=time.monotonic();ref,d,u=initialize(spec);initial=u.copy();duration=spec.get('duration_factor',1.)*ref.duration
    count=int(np.ceil(duration*ref.sound*(d.elements*(d.degree+1)**2/ref.radius+2/d.dz)/.8))*spec.get('step_factor',1);dt=duration/count
    work=0.;sources=np.zeros(6);history=[totals(d,u,ref,0.,work,sources)];snapshots=[u.copy()];snapshot_steps=[0]
    previous=u.copy();work_before=work;sources_before=sources.copy();status='complete';error=None;accepted=0
    for i in range(count):
        previous=u.copy();work_before=work;sources_before=sources.copy()
        try:result,dw,ds=step(d,u,ref,i*dt,dt,spec['manufactured'])
        except ValueError as exc:status='failed';error=repr(exc);break
        u=result;work+=dw;sources+=ds;accepted=i+1
        if accepted%max(1,count//10)==0 or accepted==count:
            history.append(totals(d,u,ref,accepted*dt,work,sources));snapshots.append(u.copy());snapshot_steps.append(accepted)
            if accepted<count:
                save(directory,f'fluid_{index:02d}',dict(spec=spec,status='running',accepted_steps=accepted,steps=count,dt=dt,source_sha256=identity['source_sha256'],elapsed_seconds=time.monotonic()-start),
                     dict(state=u,initial=initial,previous=previous,history=np.array(history),snapshots=np.array(snapshots),snapshot_steps=np.array(snapshot_steps),work_before=work_before,sources_before=sources_before))
            print(json.dumps(dict(case=index,fraction=accepted/count,elapsed_seconds=time.monotonic()-start)),flush=True)
    if snapshot_steps[-1]!=accepted:history.append(totals(d,u,ref,accepted*dt,work,sources));snapshots.append(u.copy());snapshot_steps.append(accepted)
    h=np.array(history);f=physical(u,d.r,ref);du,_=d.derivative(u)
    r,z,volume=grid(ref.radius,2*ref.half_length,32,256,24);sampled=sample(d,u,r,len(z));fi=physical(sampled,r[:,None],ref)
    integrals=np.sum(volume[...,None]*sampled,axis=(0,1));continuous=float(np.sum(volume*fi['energy']))
    metrics=(differences(sampled,ref.state(r[:,None],z[None,:],accepted*dt),r,volume,ref,ref.radius/8) if spec['manufactured'] else None)
    scales=np.sum(d.volume[...,None]*abs(initial),axis=(0,1,2));characteristic=ref.number_density*np.sum(d.volume)*np.array([1,ref.sound,ref.sound*ref.radius,ref.sound,1,1])
    drift=h[:,2:8]-h[0,2:8]-h[:,9:15];normalized=np.max(abs(drift)/np.maximum(scales,1e-12*characteristic),axis=0)
    record=dict(spec=spec,status=status,error=error,source_sha256=identity['source_sha256'],accepted_steps=accepted,steps=count,dt=dt,duration=accepted*dt,
                elapsed_seconds=time.monotonic()-start,manufactured_errors=metrics,continuous_energy=continuous,continuous_state_integrals=integrals.tolist(),
                physical_energy_change=float((h[-1,1]-h[0,1])/h[0,1]),physical_source_work=float(work/h[0,1]),
                physical_energy_work_remainder=float((h[-1,1]-h[0,1]-work)/h[0,1]),conserved_source_remainders=normalized[[0,2,3,4,5]].tolist(),
                nodal_continuous_energy_difference=float((h[-1,1]-continuous)/continuous),minimum_density=float(np.min(f['n'])),minimum_pressure=float(np.min(f['pressure'])))
    save(directory,f'fluid_{index:02d}',record,dict(initial=initial,state=u,previous=previous,rate=du,history=h,snapshots=np.array(snapshots),snapshot_steps=np.array(snapshot_steps),work_before=work_before,sources_before=sources_before))
    print(json.dumps(dict(case=index,status=status,remainder=record['physical_energy_work_remainder'],elapsed_seconds=record['elapsed_seconds'])),flush=True)
    return record


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);parser.add_argument('--workers',type=int,default=1);parser.add_argument('--cases',type=int,nargs='*');args=parser.parse_args()
    d=args.directory;identity=provenance()
    if d.exists():
        if json.loads((d/'provenance.json').read_text())['source_sha256']!=identity['source_sha256']:raise ValueError('Frozen source differs')
    else:d.mkdir(parents=True);freeze(d,identity)
    specs=protocol();indices=args.cases if args.cases is not None else range(len(specs))
    if args.workers==1:
        for i in indices:evolve(d,i,specs[i],identity)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(evolve,d,i,specs[i],identity) for i in indices]
            for future in futures:future.result()


if __name__=='__main__':main()
