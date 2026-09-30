"""Joint two-entropy original rates and unforced checkpointed axial trajectories."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import traceback
import numpy as np
from plasma_bench.two_entropy_axial import TwoEntropyAxial
from .two_entropy_reference import LogFluidReference,state,physical,parent_rate,energy_rate,measure
from .two_entropy_records import provenance,freeze,save,sha


def initialize(spec):
    kind=spec.get('kind','contact');ref=LogFluidReference(kind=kind,direction=spec.get('direction',1),gamma_i=5/3 if kind=='acoustic' else 1.4)
    solver=TwoEntropyAxial(spec['points'],2*ref.half_length,ref.mass,gammas=(ref.gamma_e,ref.gamma_i),order=spec.get('order',6))
    return ref,solver,state(ref,solver.z)


def calculate_rate(spec):
    ref,solver,_=initialize(spec);u=state(ref,solver.z,manufactured=True);actual,d=solver.derivative(u);exact,e=parent_rate(ref,solver.z,True)
    errors={key:float(np.linalg.norm(actual[:,j]-exact[:,j])/np.linalg.norm(exact[:,j]))
            for key,j in (('density',0),('radial_momentum',1),('azimuthal_momentum',2),('axial_momentum',3),('electron_entropy_density',4),('ion_entropy_density',5))}
    true_energy=energy_rate(u,actual,ref);scale=np.sum(abs(true_energy))
    result=dict(spec=spec,errors=errors,energy_rate_relative=float(np.linalg.norm(true_energy-e)/np.linalg.norm(e)),
                global_energy_rate_relative=float(np.sum(true_energy)/scale),
                local_energy_flux_relative=float(np.linalg.norm(true_energy-d['energy_flux_rate'])/np.linalg.norm(true_energy)),
                global_state_rate_relative=(np.sum(actual,axis=0)/np.sum(abs(actual),axis=0)).tolist())
    return result,dict(state=u,rate=actual,exact_rate=exact,physical_energy_rate=true_energy,exact_energy_rate=e,energy_flux_rate=d['energy_flux_rate'])


def totals(u,ref,solver,time):
    f=physical(u,ref);return np.r_[time,solver.dx*np.sum(u,axis=0),solver.dx*np.sum(f['energy']),f['n'].min(),f['pressure'].min()]


def integrate(spec,checkpoint=None):
    start=time.monotonic();ref,solver,u=initialize(spec);initial=u.copy();duration=ref.duration
    steps=4*int(np.ceil(duration*1.5*ref.sound/(.25*solver.dx)*spec.get('step_factor',1)/4));dt=duration/steps
    first,first_d=solver.derivative(u);history=[totals(u,ref,solver,0.)];times=[0.];snapshots=[u];before=u;failure=None
    def pack(status):
        record=dict(status=status,spec=spec,reference=asdict(ref),steps=steps,dt=dt,end_time=duration,
                    accepted_steps=len(history)-1,elapsed_seconds=time.monotonic()-start)
        if failure:record.update(failure)
        return record,dict(initial=initial,initial_rate=first,initial_energy_flux_rate=first_d['energy_flux_rate'],final=u,before_final=before,
                           times=np.array(times),snapshots=np.array(snapshots),history=np.array(history))
    for index in range(1,steps+1):
        previous=u
        try:u=solver.step(u,dt)
        except (ValueError,RuntimeError) as error:
            failure=dict(error=repr(error),traceback=traceback.format_exc(),attempted_step=index,last_accepted_time=(index-1)*dt);break
        before=previous;history.append(totals(u,ref,solver,index*dt))
        if index%(steps//4)==0:
            times.append(index*dt);snapshots.append(u)
            if checkpoint:checkpoint(*pack('running'))
    record,arrays=pack('failed' if failure else 'complete')
    record['norms']=[measure(v,ref,t) for v,t in zip(snapshots,times)]
    record['nodal_energy_drift']=(history[-1][7]-history[0][7])/history[0][7]
    record['continuous_energy_drift']=(record['norms'][-1]['continuous_energy']-record['norms'][0]['continuous_energy'])/record['norms'][0]['continuous_energy']
    return record,arrays


def protocol():
    return [dict(kind=kind,points=n) for kind in ('contact','acoustic') for n in (16,32,64,128)]+[
        dict(kind='contact',points=32,step_factor=2),dict(kind='contact',points=32,step_factor=4),
        dict(kind='contact',points=32,direction=-1),dict(kind='acoustic',points=32,step_factor=2),dict(kind='acoustic',points=32,direction=-1)]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);args=parser.parse_args();d=args.directory
    d.mkdir(parents=True,exist_ok=False);identity=provenance();freeze(d,identity);rates=[];trajectories=[]
    for i,spec in enumerate(dict(points=n,order=o) for o in (2,4,6) for n in (16,32,64,128)):
        r,a=calculate_rate(spec);path=d/f'rate_{i:02d}.npz';np.savez_compressed(path,**a)
        r.update(source_sha256=identity['source_sha256'],arrays_file=path.name,arrays_sha256=sha(path));rates.append(r)
        print('rate',spec,r['energy_rate_relative'],r['global_energy_rate_relative'],flush=True)
    (d/'rates.json').write_text(json.dumps(rates,indent=2)+'\n')
    for index,spec in enumerate(protocol()):
        name=f'fluid_{index:02d}'
        def persist(r,a):return save(d,name,dict(r,source_sha256=identity['source_sha256']),a)
        r,a=integrate(spec,persist);r=persist(r,a);trajectories.append(r)
        (d/'trajectories.json').write_text(json.dumps(trajectories,indent=2)+'\n')
        print(name,spec,r['status'],r['norms'][-1]['n'],'energy',r['nodal_energy_drift'],'continuous',r['continuous_energy_drift'],flush=True)
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source changed during study')


if __name__=='__main__':main()
