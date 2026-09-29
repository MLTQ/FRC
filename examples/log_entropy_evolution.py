"""Finite-time sheared log-entropy transport with physical conservation targets."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import traceback
import numpy as np
from scipy.constants import elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.variable_compatible_current import VariableCompatibleCurrent
from plasma_bench.compatible_log_thermal import CompatibleLogThermal
from plasma_bench.log_entropy_transport import LogEntropyTransport
from .log_entropy_reference import LogEntropyReference
from .finite_mass_manufactured import cell_average
from .physical_scalar_sampler import PhysicalScalarSampler
from .log_thermal_records import provenance,freeze,sha


def initialize(spec):
    reference=LogEntropyReference(toroidal=.02*spec.get('direction',1),
        amplitude=spec.get('amplitude',.3),gamma=spec.get('gamma',5/3))
    nr=spec['nr'];mesh=MHDMesh(reference.radius,reference.half_length,nr,2*nr)
    n=cell_average(mesh,lambda r,z:reference.fields(r,z)['n'],order=12)
    model=VariableCompatibleCurrent(mesh,n,reference.ion_mass,degree=spec.get('degree',5),
                                   quadrature_order=spec.get('quadrature_order',10))
    g=cell_average(mesh,lambda r,z:reference.fields(r,z)['g'],order=12)
    magnetic=model.geometry.pack(np.zeros((nr+1,2*nr+1)),g);current=model.current(magnetic)
    thermal=CompatibleLogThermal(model,gamma=reference.gamma)
    s=reference.fields(*np.meshgrid(mesh.rf,mesh.zf[:-1],indexing='ij'))['entropy']
    return reference,model,thermal,LogEntropyTransport(thermal,current),s,magnetic


def physical_norms(mesh,reference,spec,arrays,*,order=12):
    grid=PhysicalScalarSampler(mesh,spec.get('degree',5),order)
    n=grid.density(arrays['n_cell']);records=[]
    for time,state in zip(arrays['times'],arrays['snapshots']):
        s=grid.scalar(state);proxy=n*E*np.exp(s)
        pressure=proxy*(n/reference.density_scale)**(reference.gamma-1)
        internal=pressure/(reference.gamma-1)
        exact=reference.fields(grid.r,grid.z,float(time))
        temperature=pressure/(n*E);wanted_temperature=exact['pressure']/(exact['n']*E)
        background=(exact['n']/reference.density_scale)**(reference.gamma-1)*np.exp(reference.background_entropy)
        records.append(dict(time=float(time),
            entropy_wave=grid.metric(s-reference.background_entropy,exact['entropy']-reference.background_entropy),
            pressure=grid.metric(pressure,exact['pressure']),
            temperature_wave=grid.metric(temperature-background,wanted_temperature-background),
            integrals=dict(entropy=float(np.sum(grid.weight*n*s)),proxy=float(np.sum(grid.weight*proxy)),
                           internal=float(np.sum(grid.weight*internal)))))
    return records


def integrate(spec):
    ref,m,t,transport,s,magnetic=initialize(spec);geo=m.geometry
    steps=int(5*m.mesh.nz*spec.get('step_factor',1));end=ref.transit_time;dt=end/steps
    initial=s.copy();snapshots=[s.copy()];times=[0.];history=[];residual=0.;iterations=0
    def totals(value,time):
        f=t.fields(value)
        return dict(time=time,entropy=float(np.sum(geo.weight*m.density*f['entropy'])),
            proxy=float(np.sum(geo.weight*m.density*f['proxy_per_particle'])),
            internal=float(np.sum(geo.weight*f['internal'])),min_pressure=float(f['pressure'].min()),
            max_pressure=float(f['pressure'].max()))
    history.append(totals(s,0.));before_final=None
    for step in range(1,steps+1):
        if step==steps:before_final=s.copy()
        result=transport.step(s,dt);s=result.entropy
        residual=max(residual,result.max_scalar_residual);iterations+=result.iterations
        history.append(totals(s,step*dt))
        if step%(steps//4)==0:
            snapshots.append(s.copy());times.append(step*dt)
            print('transport',spec,'fraction',step/steps,flush=True)
    arrays=dict(n_cell=m.density_map.cell_means,current=transport.current,magnetic=magnetic,
                initial=initial,final=s,before_final=before_final,snapshots=np.stack(snapshots),times=np.array(times))
    norms=physical_norms(m.mesh,ref,spec,arrays,order=12)
    higher=physical_norms(m.mesh,ref,spec,arrays,order=16)
    exact=ref.integrals();conservation={}
    for key in ('entropy','proxy','internal'):
        before=history[0][key];values=np.array([h[key] for h in history])
        conservation[key]=dict(initial_relative_to_exact=float(before/exact[key]-1),
            final_relative_to_exact=float(values[-1]/exact[key]-1),
            final_relative_to_initial=float(values[-1]/before-1),
            max_relative_drift=float(np.max(abs(values/before-1))))
    record=dict(status='complete',spec=spec,reference=asdict(ref),steps=steps,dt=dt,end_time=end,
        max_scalar_residual=residual,scalar_iterations=iterations,history=history,norms=norms,
        higher_quadrature_norms=higher,exact_integrals=exact,conservation=conservation)
    return record,arrays


def protocol():
    return [dict(nr=n) for n in (8,16,32)]+[dict(nr=16,step_factor=2),dict(nr=16,quadrature_order=12),
        dict(nr=16,direction=-1),dict(nr=16,amplitude=1.),dict(nr=16,gamma=1.4),dict(nr=16,degree=3)]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path)
    directory=parser.parse_args().directory;directory.mkdir(parents=True,exist_ok=False)
    identity=provenance();freeze(directory,identity);records=[]
    for index,spec in enumerate(protocol()):
        name=f'transport_{index:02d}';print(name,spec,flush=True)
        try:
            record,arrays=integrate(spec);path=directory/(name+'.npz');np.savez_compressed(path,**arrays)
            record.update(arrays_file=path.name,arrays_sha256=sha(path))
        except Exception as error:
            record=dict(status='failed',spec=spec,error=repr(error),traceback=traceback.format_exc())
        record['source_sha256']=identity['source_sha256'];records.append(record)
        (directory/(name+'.json')).write_text(json.dumps(record,indent=2)+'\n')
        print(record['status'],record.get('norms',[{}])[-1].get('entropy_wave'),flush=True)
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Source changed during transport')
    (directory/'transports.json').write_text(json.dumps(records,indent=2)+'\n')
    if any(r['status']!='complete' for r in records):raise RuntimeError('Preserved failed transport cases')


if __name__=='__main__':main()
