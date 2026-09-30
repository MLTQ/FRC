"""Resolve the observed axis error into pressure, kinetic source and convection."""
import argparse
import json
from pathlib import Path
import numpy as np
from plasma_bench.compatible_particle_flux import CompatibleParticleFlux
from plasma_bench.compatible_bulk_thermal import CompatibleBulkThermal
from .radial_work_study import initialize,STATE_KEYS
from .radial_work_reference import RadialWorkReference
from .log_fluid_reference import initial
from .log_fluid_norms import FluidSampler
from .radial_work_records import provenance,sha


def calculate(nr):
    _,d,_=initialize(dict(kind='contact',nr=nr));ref=RadialWorkReference();state=initial(d,ref)
    rate,_=d.derivative(state);model,species,*_=d.evaluate(state);states=(state.electron,state.ion)
    flux=CompatibleParticleFlux(model.density_map);faces=flux.from_momentum(state.momentum[...,:2])
    thermal=[CompatibleBulkThermal(t,flux).rate(s,faces) for t,s in zip(species,states)]
    pressure,_=d.radial_pressure.solve(model,flux,sum(t.face_work.radial for t in thermal))
    kinetic,_=d.radial_pressure.solve(model,flux,d.radial_kinetic_work(model,flux,state))
    convection=rate.momentum[...,0]-pressure-kinetic
    grid=FluidSampler(d.mesh,5,14);n=grid.density(state.density);nt=grid.density(rate.density)
    velocity=grid.momentum(state.momentum)/n[...,None]
    def acceleration(value):
        vector=np.zeros_like(state.momentum);vector[...,0]=value
        return grid.momentum(vector)[...,0]/n
    h=1e-28;f=ref.fields(grid.r,grid.z)
    dr={k:np.imag(v)/h for k,v in ref.fields(grid.r+1j*h,grid.z).items()}
    dz={k:np.imag(v)/h for k,v in ref.fields(grid.r,grid.z+1j*h).items()}
    expected=dict(pressure=-(dr['pe']+dr['pi'])/(ref.mass*f['n']),kinetic=f['velocity'][...,1]**2/grid.r,
                  convection=-f['velocity'][...,0]*dr['velocity'][...,0]-f['velocity'][...,2]*dz['velocity'][...,0])
    actual=dict(pressure=acceleration(pressure),kinetic=acceleration(kinetic),
                convection=acceleration(convection)-velocity[...,0]*nt/n)
    metrics={k:grid.metric(v,expected[k]) for k,v in actual.items()}
    metrics['total']=grid.metric(sum(actual.values()),sum(expected.values()))
    arrays=dict(pressure_source=pressure,kinetic_source=kinetic,convection_rate=convection,density_rate=rate.density)
    arrays.update({key:getattr(state,key) for key in STATE_KEYS})
    return dict(nr=nr,metrics=metrics),arrays


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);d=parser.parse_args().directory
    d.mkdir(parents=True,exist_ok=False);records=[];identity=provenance()
    for nr in (8,16,32,64):
        record,arrays=calculate(nr);path=d/f'terms_{nr}.npz';np.savez_compressed(path,**arrays)
        record.update(arrays_file=path.name,arrays_sha256=sha(path));records.append(record)
        print(nr,{k:v['axis_relative_L2'] for k,v in record['metrics'].items()},flush=True)
    source=Path(__file__)
    (d/'records.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],probe_sha256=sha(source),records=records),indent=2)+'\n')
    (d/'probe_source.py').write_bytes(source.read_bytes());(d/'probe_source.md').write_bytes(source.with_suffix('.md').read_bytes())


if __name__=='__main__':main()
