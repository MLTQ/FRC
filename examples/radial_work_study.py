"""Unforced contact, simple-wave and rotation trajectories with preserved failures."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from dataclasses import asdict,fields
import json
from pathlib import Path
import time
import traceback
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.log_fluid import LogFluidState
from plasma_bench.radial_pressure_fluid import RadialPressureFluid
from .log_fluid_reference import LogFluidReference,initial
from .log_fluid_norms import FluidSampler
from .radial_work_records import provenance,freeze,sha

STATE_KEYS=tuple(f.name for f in fields(LogFluidState))
HISTORY_KEYS=('time','particles','axial_momentum','angular_momentum','energy_account','physical_energy',
              'electron_entropy','ion_entropy','min_density','min_electron_pressure','min_ion_pressure',
              'max_speed','spatial_energy_ledger','scalar_residual','continuity_face_agreement')


def initialize(spec):
    kind=spec['kind'];ref=LogFluidReference(kind=kind,direction=spec.get('direction',1),gamma_i=5/3 if kind=='acoustic' else 1.4)
    nr=spec['nr'];mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr)
    geometry=CompatibleGeometry(mesh,degree=spec.get('degree',5),quadrature_order=spec.get('quadrature_order',12))
    solver=RadialPressureFluid(geometry,ref.ion_mass,gamma_e=ref.gamma_e,gamma_i=ref.gamma_i)
    return ref,solver,initial(solver,ref)


def totals(solver,state,time,ledger,diagnostic=None):
    model,species,f,kinetic,energy=solver.evaluate(state);g=solver.geometry;v=solver.mesh.volumes;d=diagnostic or {}
    return [time,float(np.sum(v*state.density)),float(solver.mass*np.sum(v*state.momentum[...,1])),
            float(solver.mass*np.sum(v*state.momentum[...,2])),float(np.sum(v*state.energy_account)),float(np.sum(v*energy)),
            *[float(np.sum(g.weight*model.density*fs['entropy'])) for fs in f],float(model.density.min()),
            float(f[0]['pressure'].min()),float(f[1]['pressure'].min()),float(np.sqrt(np.sum(kinetic['velocity']**2,axis=-1)).max()),
            float(np.sum(v*ledger)),d.get('scalar_residual',0.),d.get('continuity_face_agreement',0.)]


def integrate(spec):
    start=time.monotonic();ref,solver,state=initialize(spec);mesh=solver.mesh;init=state
    end=ref.duration*spec.get('duration_factor',1.)
    steps=4*int(np.ceil(end*1.5*ref.sound/(.4*min(mesh.dr,mesh.dz))*spec.get('step_factor',1)/4));dt=end/steps
    ledger=np.zeros_like(state.density);snapshots=[state];times=[0.];ledgers=[ledger.copy()]
    physical=[solver.evaluate(state)[-1]];history=[totals(solver,state,0.,ledger)]
    first_rate,first_d=solver.derivative(state);before=state;failure=None
    for step in range(1,steps+1):
        before=state
        try:state,diagnostic=solver.step(state,(step-1)*dt,dt)
        except Exception as error:
            failure=dict(error=repr(error),traceback=traceback.format_exc(),attempted_step=step,last_accepted_time=(step-1)*dt)
            break
        ledger+=diagnostic['integrated_spatial_energy_defect']
        history.append(totals(solver,state,step*dt,ledger,diagnostic))
        if step%(steps//4)==0:
            snapshots.append(state);times.append(step*dt);ledgers.append(ledger.copy());physical.append(solver.evaluate(state)[-1])
            print('fluid',spec,'fraction',step/steps,'seconds',round(time.monotonic()-start,1),flush=True)
    arrays=dict(times=np.array(times),history=np.array(history),spatial_ledger_snapshots=np.array(ledgers),
                spatial_ledger_final=ledger,physical_energy_snapshots=np.array(physical),initial_physical_energy_rate=first_d['kinetic_rate']+sum(first_d['thermal_rates']))
    for key in STATE_KEYS:
        arrays[key+'_snapshots']=np.stack([getattr(s,key) for s in snapshots]);arrays[key+'_initial']=getattr(init,key)
        arrays[key+'_final']=getattr(state,key);arrays[key+'_before_final']=getattr(before,key);arrays[key+'_initial_rate']=getattr(first_rate,key)
    sampler=FluidSampler(mesh,spec.get('degree',5),14)
    norms=[sampler.measure(s,ref,t)[0] for s,t in zip(snapshots,times)]
    rates,independent_rate=sampler.rates(init,first_rate,ref,0.)
    agreement=float(np.linalg.norm(independent_rate-arrays['initial_physical_energy_rate'])/max(np.linalg.norm(independent_rate),1e-300))
    record=dict(status='failed' if failure else 'complete',spec=spec,reference=asdict(ref),steps=steps,dt=dt,end_time=end,
                accepted_steps=len(history)-1,history_columns=HISTORY_KEYS,norms=norms,initial_rate_norms=rates,
                independent_initial_energy_rate_relative=agreement,elapsed_seconds=time.monotonic()-start)
    if failure:record.update(failure)
    return record,arrays


def protocol():
    return [dict(kind=kind,nr=n) for kind in ('contact','acoustic','rotation') for n in (8,16,32)]+[
        dict(kind='contact',nr=16,step_factor=2),dict(kind='contact',nr=16,direction=-1),
        dict(kind='acoustic',nr=16,step_factor=2),dict(kind='acoustic',nr=16,quadrature_order=16),
        dict(kind='acoustic',nr=16,direction=-1)]


def run_case(index,spec,directory,source_sha256):
    name=f'fluid_{index:02d}';print(name,spec,flush=True)
    try:
        record,arrays=integrate(spec);path=directory/(name+'.npz');np.savez_compressed(path,**arrays)
        record.update(arrays_file=path.name,arrays_sha256=sha(path))
    except Exception as error:record=dict(status='failed',spec=spec,error=repr(error),traceback=traceback.format_exc())
    record['source_sha256']=source_sha256
    (directory/(name+'.json')).write_text(json.dumps(record,indent=2)+'\n');print(name,record['status'],flush=True)
    return index,record


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--cases',type=int,nargs='*');args=parser.parse_args();directory=args.directory
    directory.mkdir(parents=True,exist_ok=False);identity=provenance();freeze(directory,identity);completed={}
    wanted=range(len(protocol())) if args.cases is None else args.cases
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        tasks=[pool.submit(run_case,i,protocol()[i],directory,identity['source_sha256']) for i in wanted]
        for task in as_completed(tasks):
            index,record=task.result();completed[index]=record
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source changed during study')
    (directory/'trajectories.json').write_text(json.dumps([completed[i] for i in sorted(completed)],indent=2)+'\n')
    if any(r['status']!='complete' for r in completed.values()):raise RuntimeError('Preserved failed fluid trajectories')


if __name__=='__main__':main()
