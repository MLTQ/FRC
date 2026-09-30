"""Original Euler rates, separated radial forces and physical energy defects."""
import argparse
import json
from pathlib import Path
import numpy as np
from plasma_bench.compatible_particle_flux import CompatibleParticleFlux
from plasma_bench.compatible_bulk_thermal import CompatibleBulkThermal
from .face_radial_reference import initialize
from .face_radial_study import STATE_KEYS
from .face_radial_norms import FaceRadialSampler
from .radial_work_reference import RadialWorkReference
from .face_radial_records import provenance,freeze,sha


def calculate(spec):
    ref,solver,state=initialize(spec,RadialWorkReference());rate,d=solver.derivative(state)
    grid=FaceRadialSampler(solver.mesh,spec.get('degree',5),14)
    norms,energy=grid.rates(state,rate,ref,0.);physical=d['kinetic_rate']+sum(d['thermal_rates'])
    model,species,_,kinetic,initial_energy=solver.evaluate(state);g=solver.geometry;v=kinetic['velocity']
    flux=CompatibleParticleFlux(model.density_map);faces=solver.particle_faces(flux,state)
    thermal=[CompatibleBulkThermal(t,flux).rate(s,faces) for t,s in zip(species,(state.electron,state.ion))]
    pressure,_=solver.kinematics.solve_rate(model.density_map,np.zeros_like(g.r),sum(t.face_work.radial for t in thermal))
    centrifugal,_=solver.kinematics.solve_rate(model.density_map,model.density*v[...,1]**2/g.r,np.zeros_like(faces.radial))
    zero=np.zeros_like(state.momentum);f=grid.physical(state,ref);exact=ref.fields(grid.r,grid.z);dr,dz,_=ref.derivatives(grid.r,grid.z)
    actual_p=grid.momentum(pressure,zero)[...,0]/f['n'];actual_c=grid.momentum(centrifugal,zero)[...,0]/f['n']
    actual_a=(grid.momentum(rate.radial-pressure-centrifugal,zero)[...,0]-f['velocity'][...,0]*grid.density(rate.density))/f['n']
    target_p=-(dr['pe']+dr['pi'])/(ref.mass*exact['n']);target_c=exact['velocity'][...,1]**2/grid.r
    target_a=-exact['velocity'][...,0]*dr['velocity'][...,0]-exact['velocity'][...,2]*dz['velocity'][...,0]
    components={key:grid.metric(a,b) for key,a,b in (('pressure',actual_p,target_p),('centrifugal',actual_c,target_c),('advection',actual_a,target_a))}
    arrays=dict(physical_energy_rate=physical,independent_energy_rate=energy,local_energy_defect=d['local_energy_defect'],
                pressure_radial_rate=pressure,centrifugal_radial_rate=centrifugal)
    for key in STATE_KEYS:arrays[key]=getattr(state,key);arrays[key+'_rate']=getattr(rate,key)
    volumes=solver.mesh.volumes;scale=float(np.sum(volumes*initial_energy))*ref.sound/ref.radius
    return dict(spec=spec,norms=norms,radial_components=components,
                global_physical_energy_rate=float(np.sum(volumes*physical)),
                scaled_global_energy_defect=float(np.sum(volumes*d['local_energy_defect']))/scale,
                independent_energy_rate_relative=float(np.linalg.norm(energy-physical)/np.linalg.norm(energy)),
                radial_metric_residual=d['radial_metric_residual'],scalar_residual=d['scalar_residual']),arrays


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);d=parser.parse_args().directory
    d.mkdir(parents=True,exist_ok=False);identity=provenance();freeze(d,identity);records=[]
    specs=[dict(nr=n,degree=p) for p in (3,5) for n in (8,16,32,64)]+[dict(nr=32,quadrature_order=16)]
    for index,spec in enumerate(specs):
        record,arrays=calculate(spec);path=d/f'rate_{index:02d}.npz';np.savez_compressed(path,**arrays)
        record.update(arrays_file=path.name,arrays_sha256=sha(path),source_sha256=identity['source_sha256'])
        records.append(record);(d/'rates.json').write_text(json.dumps(records,indent=2)+'\n')
        print(spec,'acceleration',record['norms']['velocity']['relative_L2'],'axis',record['norms']['velocity_components']['radial']['axis_relative_L2'],
              'energy defect',record['scaled_global_energy_defect'],flush=True)
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Rate source changed')


if __name__=='__main__':main()
