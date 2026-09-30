"""Original two-dimensional rates and unforced cylindrical verification records."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.cylindrical_entropy_fluid import CylindricalEntropyFluid
from .log_fluid_reference import LogFluidReference
from .radial_work_reference import RadialWorkReference
from .cylindrical_acoustic_reference import CylindricalAcousticReference
from .cylindrical_entropy_reference import initial,physical,rates,measure
from .cylindrical_entropy_spectrum import calculate
from .cylindrical_entropy_records import provenance,freeze,save,sha


def initialize(spec):
    kind=spec.get('kind','manufactured')
    ref=(RadialWorkReference() if kind=='manufactured' else
         CylindricalAcousticReference(amplitude=spec.get('amplitude',1e-4)) if kind=='radial_acoustic' else
         LogFluidReference(kind=kind,direction=spec.get('direction',1)))
    d=CylindricalEntropyFluid(ref.radius,2*ref.half_length,spec['elements'],spec['points'],ref.mass,
                             degree=spec.get('degree',5),gammas=(ref.gamma_e,ref.gamma_i))
    return ref,d,initial(d,ref)


def rate_case(spec):
    ref,d,u=initialize(spec);du,diagnostic=d.derivative(u);actual=rates(u,du,d.r,ref);exact=ref.rates(d.r,d.z)
    targets=dict(n=exact['n'],velocity=exact['velocity'],pressure=np.stack((exact['pe'],exact['pi']),axis=-1),energy=exact['energy'])
    errors={}
    for key,wanted in targets.items():
        weight=d.volume if wanted.ndim==3 else d.volume[...,None];diff=actual[key]-wanted
        errors[key]=dict(relative_L2=float(np.sqrt(np.sum(weight*diff**2)/np.sum(weight*wanted**2))),
                         relative_max=float(np.max(abs(diff))/np.max(abs(wanted))))
    err=actual['velocity']-exact['velocity'];axis_weight=d.volume[0]
    errors['axis_radial']=dict(relative_L2=float(np.sqrt(np.sum(axis_weight*err[0,...,0]**2)/np.sum(axis_weight*exact['velocity'][0,...,0]**2))),
                               relative_max=float(np.max(abs(err[0,...,0]))/np.max(abs(exact['velocity'][0,...,0]))),
                               maximum_point_relative=float(np.max(abs(err[0,...,0]/exact['velocity'][0,...,0]))))
    sums=np.sum(d.volume[...,None]*du,axis=(0,1,2));scales=np.sum(d.volume[...,None]*abs(du),axis=(0,1,2))
    energy=actual['energy'];conservation=abs(sums[[0,2,3,4,5]])/scales[[0,2,3,4,5]]
    record=dict(spec=spec,errors=errors,conserved_rate_relative=conservation.tolist(),
                global_physical_energy_rate_relative=float(abs(np.sum(d.volume*energy))/np.sum(d.volume*abs(energy))),
                independent_energy_rate_difference=float(np.linalg.norm(energy-diagnostic['physical_energy_rate'])/np.linalg.norm(energy)),
                initial_physical_integration=measure(d,u,ref),higher_physical_integration=measure(d,u,ref,radial_order=48,axial_factor=8))
    return record,dict(state=u,rate=du,physical_energy_rate=energy,exact_velocity_rate=exact['velocity'],radius=d.r,volume=d.volume)


def protocol():
    return ([dict(kind='contact',elements=1,points=n,degree=3) for n in (16,32,64)]+
            [dict(kind='rotation',elements=n,points=8,degree=5) for n in (2,4,8)]+
            [dict(kind='radial_acoustic',elements=n,points=8,degree=3) for n in (2,4,8)]+
            [dict(kind='radial_acoustic',elements=4,points=8,degree=3,amplitude=1e-5),
             dict(kind='contact',elements=1,points=32,degree=3,step_factor=2),
             dict(kind='contact',elements=1,points=32,degree=3,direction=-1)])


def history_row(d,u,ref,t):
    f=physical(u,d.r,ref)
    return np.r_[t,np.sum(d.volume*f['energy']),np.sum(d.volume[...,None]*u,axis=(0,1,2))]


def evolve(directory,index,spec,identity):
    started=time.monotonic();ref,d,u=initialize(spec);before=u.copy()
    duration=ref.radius/ref.sound if ref.kind=='rotation' else ref.duration
    bound=ref.sound*(d.elements*(d.degree+1)**2/ref.radius+2/d.dz)
    steps=int(np.ceil(duration*bound/.8))*spec.get('step_factor',1);dt=duration/steps
    history=[history_row(d,u,ref,0.)];snapshots=[u.copy()];snapshot_steps=[0];previous=u.copy();status='complete';failure=None
    initial_metrics=measure(d,u,ref);accepted=0
    for step in range(1,steps+1):
        previous=u.copy()
        try:u=d.step(u,dt)
        except ValueError as error:status='failed';failure=str(error);break
        accepted=step
        if step%max(1,steps//10)==0 or step==steps:
            history.append(history_row(d,u,ref,step*dt));snapshots.append(u.copy());snapshot_steps.append(step)
            if step<steps:
                save(directory,f'fluid_{index:02d}',dict(spec=spec,status='running',source_sha256=identity['source_sha256'],
                     accepted_steps=step,steps=steps,dt=dt,elapsed_seconds=time.monotonic()-started),
                     dict(state=u,initial=before,previous=previous,history=np.array(history),snapshot_steps=np.array(snapshot_steps),snapshots=np.array(snapshots)))
                print(json.dumps(dict(case=index,fraction=step/steps,elapsed_seconds=time.monotonic()-started)),flush=True)
    if snapshot_steps[-1]!=accepted:
        history.append(history_row(d,u,ref,accepted*dt));snapshots.append(u.copy());snapshot_steps.append(accepted)
    du,diag=d.derivative(u);final=measure(d,u,ref,accepted*dt);high=measure(d,u,ref,accepted*dt,48,8)
    f=physical(u,d.r,ref);f0=physical(before,d.r,ref);h=np.array(history)
    scale=np.sum(d.volume[...,None]*abs(before),axis=(0,1,2));fallback=ref.number_density*np.sum(d.volume)*np.array([1,ref.sound,ref.radius*ref.sound,ref.sound,1,1])
    record=dict(spec=spec,status=status,failure=failure,source_sha256=identity['source_sha256'],accepted_steps=accepted,steps=steps,dt=dt,
                elapsed_seconds=time.monotonic()-started,duration=accepted*dt,initial_metrics=initial_metrics,final_metrics=final,higher_metrics=high,
                physical_energy_drift=float((h[-1,1]-h[0,1])/h[0,1]),
                history_conserved_drift=np.max(abs(h[:,2:]-h[0,2:])/np.maximum(scale,1e-12*fallback),axis=0)[[0,2,3,4,5]].tolist(),
                max_specific_entropy_change=float(np.max(abs(f['entropy']-f0['entropy']))),
                maximum_radial_velocity_over_sound=float(np.max(abs(f['velocity'][...,0]))/ref.sound))
    if ref.kind=='radial_acoustic':
        exact=ref.fields(d.r,d.z,accepted*dt)
        record['density_perturbation_error']=float(np.sqrt(np.sum(d.volume*(f['n']-exact['n'])**2)/np.sum(d.volume*(before[...,0]-ref.number_density)**2)))
        record['radial_velocity_error_over_amplitude']=float(np.sqrt(np.sum(d.volume*(f['velocity'][...,0]-exact['velocity'][...,0])**2)/np.sum(d.volume))/(ref.amplitude*ref.sound))
    save(directory,f'fluid_{index:02d}',record,dict(state=u,initial=before,previous=previous,rate=du,history=h,snapshot_steps=np.array(snapshot_steps),snapshots=np.array(snapshots)))
    print(json.dumps(dict(case=index,status=status,energy_drift=record['physical_energy_drift'],elapsed_seconds=record['elapsed_seconds'])),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);parser.add_argument('--mode',choices=('rates','spectra','evolution'),required=True);parser.add_argument('--cases',type=int,nargs='*')
    args=parser.parse_args();d=args.directory;identity=provenance()
    if d.exists():
        if json.loads((d/'provenance.json').read_text())['source_sha256']!=identity['source_sha256']:raise ValueError('Frozen study source changed')
    else:d.mkdir(parents=True);freeze(d,identity)
    if args.mode=='evolution':
        specs=protocol()
        for i in args.cases if args.cases is not None else range(len(specs)):
            if (d/f'fluid_{i:02d}.json').exists():raise ValueError('Terminal record already exists')
            evolve(d,i,specs[i],identity)
    else:
        specs=([dict(elements=n,points=128,degree=p) for p in (3,5,7) for n in (2,4,8,16,32)] if args.mode=='rates' else
               [dict(elements=n,degree=p) for p in (3,5) for n in (2,4,8)])
        for i,spec in enumerate(specs):
            stem=f'{args.mode}_{i:02d}';path=d/(stem+'.npz')
            if path.exists():raise ValueError('Terminal record already exists')
            record,arrays=(rate_case(spec) if args.mode=='rates' else calculate(spec));np.savez_compressed(path,**arrays)
            record.update(source_sha256=identity['source_sha256'],arrays_file=path.name,arrays_sha256=sha(path))
            (d/(stem+'.json')).write_text(json.dumps(record,indent=2)+'\n');print(stem,flush=True)


if __name__=='__main__':main()
