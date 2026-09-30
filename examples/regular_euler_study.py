"""Original regular-coordinate rates and matched mixed-flow trajectories."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from .regular_euler_reference import RegularEulerReference,initial,energy
from .regular_euler_sampling import sample
from .regular_euler_step import step,integrals
from .regular_euler_spectrum import calculate
from .regular_euler_records import provenance,freeze,save,sha
from .mixed_euler_study import protocol
from .mixed_euler_sampling import grid,differences
from .cylindrical_entropy_reference import physical


def initialize(spec):
    ref=RegularEulerReference();d=RegularCoordinateFluid(ref.radius,2*ref.half_length,spec['elements'],spec['points'],ref.mass,degree=spec.get('degree',5))
    return ref,d,initial(d,ref)


def totals(d,u,ref,t,work,sources):return np.r_[t,np.sum(d.volume*energy(u,d.x,ref)),integrals(d,u),work,sources]


def rate_case(spec):
    ref,d,u=initialize(spec);du,diag=d.derivative(u);wanted=ref.parent_rate(d.x,d.z);wanted[0,0,:,1:3]=0.;diff=du-wanted
    errors=np.sqrt(np.sum(d.volume[...,None]*diff*diff,axis=(0,1,2))/np.sum(d.volume[...,None]*wanted*wanted,axis=(0,1,2)))
    axis=float(np.sqrt(np.sum(d.volume[0,1:]*diff[0,1:,:,1]**2)/np.sum(d.volume[0,1:]*wanted[0,1:,:,1]**2)))
    et=np.imag(energy(u+1j*1e-28*du,d.x,ref))/1e-28;angular=2*d.x*du[...,2]
    record=dict(spec=spec,original_state_rate_errors=errors.tolist(),first_element_active_q_error=axis,
                physical_energy_rate_defect=float(abs(np.sum(d.volume*et))/np.sum(d.volume*abs(et))),
                local_energy_chain_rule_difference=float(np.linalg.norm(et-diag['physical_energy_rate'])/np.linalg.norm(et)),
                angular_rate_defect=float(abs(np.sum(d.volume*angular))/np.sum(d.volume*abs(angular))))
    return record,dict(state=u,rate=du,exact_rate=wanted,physical_energy_rate=et,volume=d.volume,x=d.x)


def evolve(directory,index,spec,identity):
    directory=Path(directory)
    if (directory/f'fluid_{index:02d}.json').exists():raise ValueError('Terminal result already exists')
    start=time.monotonic();ref,d,u=initialize(spec);first=u.copy();duration=spec.get('duration_factor',1.)*ref.duration
    count=int(np.ceil(duration/ref.duration*(1.5*d.elements*(d.degree+1)**2+2*ref.radius/d.dz)/.8))*spec.get('step_factor',1);dt=duration/count
    work=0.;sources=np.zeros(6);history=[totals(d,u,ref,0.,work,sources)];snapshots=[u.copy()];snapshot_steps=[0]
    previous=u.copy();work_before=0.;sources_before=sources.copy();status='complete';error=None;accepted=0
    for i in range(count):
        previous=u.copy();work_before=work;sources_before=sources.copy()
        try:result,dw,ds=step(d,u,ref,i*dt,dt,spec['manufactured'])
        except ValueError as exc:status='failed';error=repr(exc);break
        u=result;work+=dw;sources+=ds;accepted=i+1
        if accepted%max(1,count//10)==0 or accepted==count:
            history.append(totals(d,u,ref,accepted*dt,work,sources));snapshots.append(u.copy());snapshot_steps.append(accepted)
            if accepted<count:save(directory,f'fluid_{index:02d}',dict(spec=spec,status='running',accepted_steps=accepted,steps=count,dt=dt,source_sha256=identity['source_sha256'],elapsed_seconds=time.monotonic()-start),
                dict(state=u,initial=first,previous=previous,history=np.array(history),snapshots=np.array(snapshots),snapshot_steps=np.array(snapshot_steps),work_before=work_before,sources_before=sources_before))
            print(json.dumps(dict(case=index,fraction=accepted/count,elapsed_seconds=time.monotonic()-start)),flush=True)
    if snapshot_steps[-1]!=accepted:history.append(totals(d,u,ref,accepted*dt,work,sources));snapshots.append(u.copy());snapshot_steps.append(accepted)
    h=np.array(history);f=d.fields(u);du,_=d.derivative(u)
    r,z,volume=grid(ref.radius,2*ref.half_length,32,256,24);sampled=sample(d,u,r,len(z));fi=physical(sampled,r[:,None],ref);continuous=float(np.sum(volume*fi['energy']))
    metrics=differences(sampled,ref.state(r[:,None],z[None,:],accepted*dt),r,volume,ref,ref.radius/8) if spec['manufactured'] else None
    scales=integrals(d,abs(first));characteristic=ref.number_density*np.sum(d.volume)*np.array([1,ref.sound/ref.radius,ref.sound*ref.radius,ref.sound,1,1])
    drift=h[:,2:8]-h[0,2:8]-h[:,9:15];normalized=np.max(abs(drift)/np.maximum(scales,1e-12*characteristic),axis=0)
    record=dict(spec=spec,status=status,error=error,source_sha256=identity['source_sha256'],accepted_steps=accepted,steps=count,dt=dt,duration=accepted*dt,elapsed_seconds=time.monotonic()-start,
        manufactured_errors=metrics,continuous_energy=continuous,physical_energy_change=float((h[-1,1]-h[0,1])/h[0,1]),physical_source_work=float(work/h[0,1]),
        physical_energy_work_remainder=float((h[-1,1]-h[0,1]-work)/h[0,1]),conserved_source_remainders=normalized[[0,2,3,4,5]].tolist(),
        nodal_continuous_energy_difference=float((h[-1,1]-continuous)/continuous),maximum_sampled_velocity_over_sound=float(np.max(abs(fi['velocity']))/ref.sound),
        minimum_density=float(np.min(f['n'])),minimum_pressure=float(np.min(f['pressure'])))
    save(directory,f'fluid_{index:02d}',record,dict(initial=first,state=u,previous=previous,rate=du,history=h,snapshots=np.array(snapshots),snapshot_steps=np.array(snapshot_steps),work_before=work_before,sources_before=sources_before))
    print(json.dumps(dict(case=index,status=status,remainder=record['physical_energy_work_remainder'],elapsed_seconds=record['elapsed_seconds'])),flush=True);return record


def main():
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--mode',choices=('rates','spectra','evolution'),required=True);p.add_argument('--workers',type=int,default=1);p.add_argument('--cases',type=int,nargs='*');args=p.parse_args()
    d=args.directory;identity=provenance()
    if d.exists():
        if json.loads((d/'provenance.json').read_text())['source_sha256']!=identity['source_sha256']:raise ValueError('Frozen source mismatch')
    else:d.mkdir(parents=True);freeze(d,identity)
    if args.mode=='evolution':
        specs=protocol();indices=args.cases if args.cases is not None else range(len(specs))
        if args.workers==1:
            for i in indices:evolve(d,i,specs[i],identity)
        else:
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                futures=[pool.submit(evolve,d,i,specs[i],identity) for i in indices]
                for future in futures:future.result()
    else:
        specs=[dict(elements=n,points=128,degree=p) for p in (3,5) for n in (2,4,8,16)]
        for i,spec in enumerate(specs):
            path=d/f'{args.mode}_{i:02d}.npz'
            if path.exists():raise ValueError('Terminal evidence already exists')
            record,arrays=rate_case(spec) if args.mode=='rates' else calculate(spec['elements'],spec['degree'])
            np.savez_compressed(path,**arrays);record.update(source_sha256=identity['source_sha256'],arrays_file=path.name,arrays_sha256=sha(path))
            path.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n');print(path.name,flush=True)


if __name__=='__main__':main()
