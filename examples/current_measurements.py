"""Reproducible complete current-inertia component measurements and archives."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
from time import perf_counter
import numpy as np
import scipy
from scipy.constants import mu_0 as MU0
from plasma_bench.moving_field import MovingField
from .current_cases import protocol
from .current_reference_integrals import integrate
from .verification_records import write_record


def provenance():
    files=sorted(Path('plasma_bench').glob('*.py'))+[Path(f'examples/{name}.py') for name in
           ('current_cases','current_measurements','current_reference_integrals','inertial_sphere_reference',
            'inertial_screening_reference','profiled_current_reference','spherical_decay_reference','verification_records')]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    return {'source_sha256':hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest(),
            'source_files':hashes,'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()}


def snapshot(o,a,bias,heat):
    d=o.domain;u=o.induced(a,bias);q=o.current_over_radius(a,bias)
    e=o.elements[0];nc=d.fluid_cell_count;br,bz=MovingField.fields(e,u,bias)
    return {'coordinates':d.fluid_reference.copy(),'cells':d.cells[:nc].copy(),'material_a':a.copy(),
            'induced_a':u,'current_over_radius':q,'current_quadrature':o.current_quadrature(q),
            'Br':br[:nc],'Bz':bz[:nc],'quadrature_coordinates':e.points[:nc].copy(),
            'volume_weights':2*np.pi*e.weights[:nc]*e.points[:nc,:,0],
            'electron_length_squared':o.de_squared.copy(),'resistivity':o.resistivity.copy(),
            'cell_heat':heat.copy(),'bias':np.asarray(bias)}


def errors(ref,kind,data,time):
    r,z=data['quadrature_coordinates'][...,0],data['quadrature_coordinates'][...,1]
    w=data['volume_weights'];exact=ref.fields(r,z) if kind=='screening' else ref.fields(r,z,time)
    def relative(value,target):return float(np.sqrt(np.sum(w*(value-target)**2)/np.sum(w*target**2)))
    berror=np.sum(w*((data['Br']-exact['Br'])**2+(data['Bz']-exact['Bz'])**2))/np.sum(w*(exact['Br']**2+exact['Bz']**2))
    result={'field_relative_volume_L2':float(np.sqrt(berror)),
            'current_relative_volume_L2':relative(data['current_quadrature'],exact['Jtheta'])}
    if kind!='screening' and time>0:
        target=np.sum(w*ref.heat_density(r,z,time),axis=1);volume=w.sum(axis=1)
        result['local_heat_relative_volume_L2']=float(np.sqrt(np.sum((data['cell_heat']-target)**2/volume)/np.sum(target**2/volume)))
    return result


def measure(case,directory):
    directory=Path(directory);source=provenance();start=perf_counter();o,ref,a=case.build()
    heat=np.zeros(o.domain.fluid_cell_count);work=0.;nsteps=int(np.ceil(case.end_time/case.dt))
    times=np.linspace(0.,case.end_time,nsteps+1);initial=snapshot(o,a,0.,heat)
    eb0,ec0=o.energies(initial['induced_a'],MU0*initial['current_over_radius'])
    arrays={f'initial_{k}':v for k,v in initial.items()};history=[]
    for left,right in zip(times[:-1],times[1:]):
        b0=ref.bias*left/case.end_time if case.kind=='screening' else 0.
        b1=ref.bias*right/case.end_time if case.kind=='screening' else 0.
        step=o.step(a,right-left,bias_start=b0,bias_end=b1);a=step.material_a;heat+=step.cell_heat;work+=step.source_work
        history.append({'time_s':float(right),'magnetic_energy_j':step.final_magnetic_energy,'current_kinetic_energy_j':step.final_current_energy,
                        'ohmic_heat_j':float(heat.sum()),'source_work_j':work,
                        'response_residual_j':float(step.final_magnetic_energy+step.final_current_energy+heat.sum()-work-eb0-ec0),
                        'step_balance_residual_j':step.balance_residual,'ampere_relative_residual':step.ampere_relative_residual,
                        'ohm_relative_residual':step.ohm_relative_residual,'vacuum_relative_residual':step.vacuum_relative_residual})
    final=snapshot(o,a,b1,heat);arrays.update({f'final_{k}':v for k,v in final.items()})
    low=integrate(ref,case.kind,case.end_time,48);high=integrate(ref,case.kind,case.end_time,96)
    continuum={**high,'quadrature_relative_change':max(abs(low[k]-high[k])/max(abs(high[k]),np.finfo(float).tiny) for k in high)}
    if case.kind=='decay':
        exact=ref.integrals(case.end_time)
        for measured,analytic in (('magnetic_energy_j','induced_energy_j'),('current_kinetic_energy_j','current_kinetic_energy_j'),('ohmic_heat_j','ohmic_heat_j')):
            if abs(high[measured]-exact[analytic])>1e-11*exact['total_response_energy_j']:
                raise ValueError('Independent sphere volume integral differs from closed energy')
    row={'name':case.name,'completed':True,'configuration':case.configuration(),'provenance':source,
         'elapsed_seconds':perf_counter()-start,'initial':{'magnetic_energy_j':eb0,'current_kinetic_energy_j':ec0},
         'initial_errors':None if case.kind=='screening' else errors(ref,case.kind,initial,0.),
         'final_errors':errors(ref,case.kind,final,case.end_time),'continuum_final':continuum,'history':history,'final':history[-1]}
    if provenance()!=source:raise RuntimeError('Source changed during current measurement')
    directory.mkdir(parents=True,exist_ok=True);path=directory/f'{case.name}.npz';np.savez_compressed(path,**arrays)
    row['fields']={'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                   'units':'SI: a T, coordinates m, B T, current A/m², q A/m³, heat J, volume m³, de² m², resistivity ohm m'}
    write_record(directory/f'{case.name}.json',row);print('Completed',case.name,flush=True)
    return row


def run(directory='outputs/current_verification',reuse=False):
    directory=Path(directory);source=provenance();rows=[]
    for case in protocol():
        path=directory/f'{case.name}.json';row=None
        if reuse and path.exists():
            old=json.loads(path.read_text());fields=directory/old['fields']['file']
            if old.get('completed') and old['configuration']==case.configuration() and old['provenance']==source and hashlib.sha256(fields.read_bytes()).hexdigest()==old['fields']['sha256']:
                row=old;print('Reuse',case.name,flush=True)
        rows.append(measure(case,directory) if row is None else row)
    write_record(directory/'measurements.json',{'completed':True,'provenance':source,'cases':rows})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--reuse',action='store_true');args=parser.parse_args();run(reuse=args.reuse)
