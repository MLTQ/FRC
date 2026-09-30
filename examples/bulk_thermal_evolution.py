"""Changing-density compression and full-transit two-temperature advection."""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from dataclasses import asdict
import json
from pathlib import Path
import traceback
import numpy as np
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.compatible_geometry import CompatibleGeometry
from plasma_bench.bulk_entropy_transport import BulkEntropyState,PrescribedBulkTransport
from .bulk_compression_reference import BulkCompressionReference
from .finite_mass_manufactured import cell_average
from .bulk_thermal_norms import transport_norms
from .bulk_thermal_records import provenance,freeze,sha


class BulkEvolutionFailure(RuntimeError):
    def __init__(self,error,record,arrays):
        super().__init__(str(error));self.record=record;self.arrays=arrays


def initialize(spec):
    ref=BulkCompressionReference(kind=spec.get('kind','compression'),speed=2000*spec.get('direction',1))
    nr=spec['nr'];mesh=MHDMesh(ref.radius,ref.half_length,nr,2*nr)
    geo=CompatibleGeometry(mesh,degree=spec.get('degree',5),quadrature_order=spec.get('quadrature_order',12))
    n=cell_average(mesh,lambda r,z:ref.fields(r,z)['n'],order=16)
    r,z=np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij');f=ref.fields(r,z)
    state=BulkEntropyState(n,f['electron'],f['ion'])
    solver=PrescribedBulkTransport(geo,ref.ion_mass,ref.momentum,gamma_e=ref.gamma_e,gamma_i=ref.gamma_i)
    return ref,solver,state


def totals(solver,state,time):
    model,species,fields=solver.evaluate(state);geo=model.geometry
    result=dict(time=float(time),particles=float(np.sum(model.mesh.volumes*state.density)),min_density=float(model.density.min()))
    for label,f in zip(('electron','ion'),fields):
        result[label]=dict(entropy=float(np.sum(geo.weight*model.density*f['entropy'])),internal=float(np.sum(geo.weight*f['internal'])),
                           min_pressure=float(f['pressure'].min()),max_pressure=float(f['pressure'].max()))
    return result


def integrate(spec):
    ref,solver,state=initialize(spec);mesh=solver.mesh;initial=state
    steps=int((2 if ref.kind=='compression' else 5)*mesh.nz*spec.get('step_factor',1));dt=ref.duration/steps
    snapshots=[state];times=[0.];history=[totals(solver,state,0.)];work=np.zeros(2);work_history=[work.copy()];residual=0.
    for step in range(1,steps+1):
        before=state
        try:result=solver.step(state,(step-1)*dt,dt)
        except Exception as error:
            retained=dict(times=np.array(times),work_history=np.array(work_history))
            for key in ('density','electron','ion'):
                retained[key+'_snapshots']=np.stack([getattr(s,key) for s in snapshots])
                retained[key+'_initial']=getattr(initial,key);retained[key+'_last_accepted']=getattr(state,key)
            detail=dict(status='failed',spec=spec,reference=asdict(ref),attempted_step=step,dt=dt,
                        last_accepted_time=(step-1)*dt,history=history,original_error=repr(error))
            raise BulkEvolutionFailure(error,detail,retained) from error
        state=result.state;work+=result.thermal_work
        residual=max(residual,result.max_scalar_residual)
        history.append(totals(solver,state,step*dt));work_history.append(work.copy())
        if step%(steps//4)==0:
            snapshots.append(state);times.append(step*dt);print('transport',spec,'fraction',step/steps,flush=True)
    arrays=dict(times=np.array(times),work_history=np.array(work_history))
    for key in ('density','electron','ion'):
        arrays[key+'_snapshots']=np.stack([getattr(s,key) for s in snapshots])
        arrays[key+'_initial']=getattr(initial,key);arrays[key+'_final']=getattr(state,key);arrays[key+'_before_final']=getattr(before,key)
    norms=transport_norms(mesh,ref,spec,arrays,order=14);higher=transport_norms(mesh,ref,spec,arrays,order=16)
    conservation=dict(particles=max(abs(h['particles']/history[0]['particles']-1) for h in history))
    for index,label in enumerate(('electron','ion')):
        initial_h=history[0][label]['internal'];values=np.array([h[label]['internal'] for h in history])
        mismatch=values-initial_h-arrays['work_history'][:,index]
        conservation[label]=dict(entropy_max_relative_drift=max(abs(h[label]['entropy']/history[0][label]['entropy']-1) for h in history),
            work_max_relative_defect=float(np.max(abs(mismatch))/initial_h),work_final_relative_defect=float(mismatch[-1]/initial_h),
            physical_energy_change=float(values[-1]-initial_h),integrated_face_work=float(work[index]))
    return dict(status='complete',spec=spec,reference=asdict(ref),steps=steps,dt=dt,end_time=ref.duration,history=history,
        max_scalar_residual=residual,norms=norms,higher_quadrature_norms=higher,conservation=conservation),arrays


def protocol():
    return [dict(nr=n,kind=kind) for kind in ('compression','translation') for n in (8,16,32)]+[
        dict(nr=16,kind='compression',step_factor=2),dict(nr=16,kind='compression',quadrature_order=16),
        dict(nr=16,kind='compression',direction=-1),dict(nr=16,kind='compression',degree=3)]


def run_case(index,spec,directory,source_sha256):
    name=f'transport_{index:02d}';print(name,spec,flush=True)
    try:
        record,arrays=integrate(spec);path=directory/(name+'.npz');np.savez_compressed(path,**arrays)
        record.update(arrays_file=path.name,arrays_sha256=sha(path))
    except Exception as error:
        record=dict(status='failed',spec=spec,error=repr(error),traceback=traceback.format_exc())
        if isinstance(error,BulkEvolutionFailure):
            record.update(error.record);path=directory/(name+'_failed_state.npz');np.savez_compressed(path,**error.arrays)
            record.update(arrays_file=path.name,arrays_sha256=sha(path))
    record['source_sha256']=source_sha256
    (directory/(name+'.json')).write_text(json.dumps(record,indent=2)+'\n');print(name,record['status'],flush=True)
    return index,record


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);parser.add_argument('--workers',type=int,choices=(1,2),default=2)
    args=parser.parse_args();directory=args.directory
    directory.mkdir(parents=True,exist_ok=False);identity=provenance();freeze(directory,identity)
    completed={}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_case,i,spec,directory,identity['source_sha256']) for i,spec in enumerate(protocol())]
        for future in as_completed(futures):
            index,record=future.result();completed[index]=record
    records=[completed[i] for i in range(len(protocol()))]
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source changed during evolution')
    (directory/'transports.json').write_text(json.dumps(records,indent=2)+'\n')
    if any(r['status']!='complete' for r in records):raise RuntimeError('Preserved failed evolution cases')


if __name__=='__main__':main()
