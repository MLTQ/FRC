"""Original two-dimensional Euler rate convergence for the radial work variant."""
import argparse
import json
from pathlib import Path
import numpy as np
from .radial_work_study import initialize,STATE_KEYS
from .log_fluid_reference import initial
from .radial_work_reference import RadialWorkReference
from .log_fluid_norms import FluidSampler
from .radial_work_records import provenance,freeze,sha


def calculate(spec):
    _,solver,_=initialize(dict(kind='contact',**spec));ref=RadialWorkReference();state=initial(solver,ref)
    rate,d=solver.derivative(state);grid=FluidSampler(solver.mesh,spec.get('degree',5),14)
    norms,energy=grid.rates(state,rate,ref,0.);physical=d['kinetic_rate']+sum(d['thermal_rates'])
    arrays=dict(physical_energy_rate=physical,independent_energy_rate=energy,local_energy_defect=d['local_energy_defect'])
    for key in STATE_KEYS:arrays[key]=getattr(state,key);arrays[key+'_rate']=getattr(rate,key)
    return dict(spec=spec,norms=norms,independent_energy_rate_relative=float(np.linalg.norm(energy-physical)/np.linalg.norm(energy)),
                radial_pressure_residual=d['radial_pressure_residual'],scalar_residual=d['scalar_residual']),arrays


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);d=parser.parse_args().directory
    d.mkdir(parents=True,exist_ok=False);identity=provenance();freeze(d,identity);records=[]
    specs=[dict(nr=n,degree=p) for p in (3,5) for n in (8,16,32,64)]+[dict(nr=32,quadrature_order=16)]
    for index,spec in enumerate(specs):
        record,arrays=calculate(spec);path=d/f'rate_{index:02d}.npz';np.savez_compressed(path,**arrays)
        record.update(arrays_file=path.name,arrays_sha256=sha(path),source_sha256=identity['source_sha256'])
        records.append(record);print(spec,record['norms']['velocity'],flush=True)
    if provenance()['source_sha256']!=identity['source_sha256']:raise RuntimeError('Rate source changed')
    (d/'rates.json').write_text(json.dumps(records,indent=2)+'\n')


if __name__=='__main__':main()
