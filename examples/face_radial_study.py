"""Checkpointed unforced trajectories for the physical common-radial-field fluid."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from dataclasses import asdict,fields
import json
from pathlib import Path
import time
import traceback
import numpy as np
from plasma_bench.face_radial_fluid import FaceRadialState
from .face_radial_reference import initialize
from .face_radial_norms import FaceRadialSampler
from .face_radial_records import provenance,freeze,sha,save
from .radial_work_study import HISTORY_KEYS

STATE_KEYS=tuple(f.name for f in fields(FaceRadialState))


def totals(solver,state,time,ledger,diagnostic=None):
    model,_,f,kinetic,energy=solver.evaluate(state);g=solver.geometry;v=solver.mesh.volumes;d=diagnostic or {}
    return [time,float(np.sum(v*state.density)),float(solver.mass*np.sum(v*state.momentum[...,0])),
            float(solver.mass*np.sum(v*state.momentum[...,1])),float(np.sum(v*state.energy_account)),float(np.sum(v*energy)),
            *[float(np.sum(g.weight*model.density*fs['entropy'])) for fs in f],float(model.density.min()),
            float(f[0]['pressure'].min()),float(f[1]['pressure'].min()),float(np.sqrt(np.sum(kinetic['velocity']**2,axis=-1)).max()),
            float(np.sum(v*ledger)),d.get('scalar_residual',0.),d.get('continuity_face_agreement',0.)]


def integrate(spec,checkpoint=None,resume=None):
    start=time.monotonic();ref,solver,init=initialize(spec);mesh=solver.mesh;state=init
    end=ref.duration*spec.get('duration_factor',1.)
    steps=4*int(np.ceil(end*1.5*ref.sound/(.4*min(mesh.dr,mesh.dz))*spec.get('step_factor',1)/4));dt=end/steps
    ledger=np.zeros_like(state.density);snapshots=[state];times=[0.];ledgers=[ledger.copy()]
    physical=[solver.evaluate(state)[-1]];history=[totals(solver,state,0.,ledger)]
    first_rate,first_d=solver.derivative(init);before=state;failure=None;elapsed=0.
    if resume is not None:
        record,a=resume
        if record['spec']!=spec or record['steps']!=steps or record['dt']!=dt:raise RuntimeError('Checkpoint input mismatch')
        def load(suffix,index=None):return FaceRadialState(*(a[k+suffix] if index is None else a[k+suffix][index] for k in STATE_KEYS))
        if not all(np.array_equal(getattr(init,k),a[k+'_initial']) for k in STATE_KEYS):raise RuntimeError('Checkpoint platform/initial-state mismatch')
        if not all(np.array_equal(getattr(first_rate,k),a[k+'_initial_rate']) for k in STATE_KEYS):raise RuntimeError('Checkpoint platform/rate mismatch')
        state=load('_final');before=load('_before_final');ledger=a['spatial_ledger_final'].copy()
        history=list(a['history']);times=list(a['times']);snapshots=[load('_snapshots',i) for i in range(len(times))]
        ledgers=list(a['spatial_ledger_snapshots']);physical=list(a['physical_energy_snapshots']);elapsed=record['elapsed_seconds']

    def bundle(status):
        arrays=dict(times=np.array(times),history=np.array(history),spatial_ledger_snapshots=np.array(ledgers),
                    spatial_ledger_final=ledger,physical_energy_snapshots=np.array(physical),
                    initial_physical_energy_rate=first_d['kinetic_rate']+sum(first_d['thermal_rates']))
        for key in STATE_KEYS:
            arrays[key+'_snapshots']=np.stack([getattr(s,key) for s in snapshots]);arrays[key+'_initial']=getattr(init,key)
            arrays[key+'_final']=getattr(state,key);arrays[key+'_before_final']=getattr(before,key);arrays[key+'_initial_rate']=getattr(first_rate,key)
        record=dict(status=status,spec=spec,reference=asdict(ref),steps=steps,dt=dt,end_time=end,
                    accepted_steps=len(history)-1,history_columns=HISTORY_KEYS,elapsed_seconds=elapsed+time.monotonic()-start)
        if failure:record.update(failure)
        return record,arrays

    for step in range(len(history),steps+1):
        previous=state
        try:state,diagnostic=solver.step(state,(step-1)*dt,dt)
        except Exception as error:
            failure=dict(error=repr(error),traceback=traceback.format_exc(),attempted_step=step,last_accepted_time=(step-1)*dt)
            break
        before=previous;ledger+=diagnostic['integrated_spatial_energy_defect']
        history.append(totals(solver,state,step*dt,ledger,diagnostic))
        if step%(steps//4)==0:
            snapshots.append(state);times.append(step*dt);ledgers.append(ledger.copy());physical.append(solver.evaluate(state)[-1])
            print('face fluid',spec,'fraction',step/steps,'seconds',round(elapsed+time.monotonic()-start,1),flush=True)
        if checkpoint and (step%max(1,steps//20)==0 or step==steps):checkpoint(*bundle('running'))
    record,arrays=bundle('failed' if failure else 'complete')
    sampler=FaceRadialSampler(mesh,spec.get('degree',5),14)
    record['norms']=[sampler.measure(s,ref,t)[0] for s,t in zip(snapshots,times)]
    record['initial_rate_norms'],independent=sampler.rates(init,first_rate,ref,0.)
    record['independent_initial_energy_rate_relative']=float(np.linalg.norm(independent-arrays['initial_physical_energy_rate'])/max(np.linalg.norm(independent),1e-300))
    return record,arrays


def protocol():
    return [dict(kind=kind,nr=n,step_factor=2 if kind=='contact' else 1)
            for kind in ('contact','acoustic','rotation') for n in (8,16,32)]+[
        dict(kind='contact',nr=8,step_factor=4),dict(kind='contact',nr=8,step_factor=2,direction=-1),
        dict(kind='acoustic',nr=16,step_factor=2),dict(kind='acoustic',nr=16,quadrature_order=16),
        dict(kind='acoustic',nr=16,direction=-1)]


def run_case(index,spec,directory,identity,resume_allowed=False):
    name=f'fluid_{index:02d}';terminal=directory/(name+'.json');checkpoint=directory/(name+'_checkpoint.json')
    if terminal.exists():return index,json.loads(terminal.read_text())
    resume=None
    if resume_allowed and checkpoint.exists():
        r=json.loads(checkpoint.read_text());path=directory/r['arrays_file']
        if r['source_sha256']!=identity or sha(path)!=r['arrays_sha256']:raise RuntimeError('Checkpoint hash/source mismatch')
        with np.load(path) as a:resume=(r,{k:a[k].copy() for k in a.files})
    def persist(record,arrays):return save(directory,name,dict(record,source_sha256=identity),arrays)
    print(name,spec,'resume',resume[0]['accepted_steps'] if resume else 0,flush=True)
    record,arrays=integrate(spec,persist,resume);record=persist(record,arrays)
    print(name,record['status'],flush=True);return index,record


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--cases',type=int,nargs='*');parser.add_argument('--resume',action='store_true');args=parser.parse_args();d=args.directory
    identity=provenance()
    if args.resume:
        if json.loads((d/'provenance.json').read_text())['source_sha256']!=identity['source_sha256']:raise RuntimeError('Resume source mismatch')
    else:d.mkdir(parents=True,exist_ok=False);freeze(d,identity)
    wanted=range(len(protocol())) if args.cases is None else args.cases;completed={}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        tasks=[pool.submit(run_case,i,protocol()[i],d,identity['source_sha256'],args.resume) for i in wanted]
        for task in as_completed(tasks):
            index,record=task.result();completed[index]=record
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source changed during study')
    (d/'trajectories.json').write_text(json.dumps([completed[i] for i in sorted(completed)],indent=2)+'\n')
    if any(r['status']!='complete' for r in completed.values()):raise RuntimeError('Preserved failed fluid trajectories')


if __name__=='__main__':main()
