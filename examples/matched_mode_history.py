"""Independent reference-parallel and orthogonal physical response histories."""
import json
from pathlib import Path
import numpy as np
from .matched_mode_comparison import load
from .matched_mode_fixture import fixture
from .matched_mode_physical import signal
from .matched_mode_rates import fingerprint
from .current_transport_sampling import TransportSampler
from .implicit_plasma_records import provenance,freeze,sha


NAMES=['density','radial_velocity','azimuthal_velocity','axial_velocity','electron_entropy','ion_entropy','radial_B','azimuthal_B','axial_B','electron_pressure','ion_pressure','current']


def calculate(audit_directory,orders=(12,16)):
    audit,_=load(audit_directory);directories=[Path(p) for p in audit['input_directories']]
    records=[json.loads((p/'summary.json').read_text()) for p in directories]
    mode,d,_,system,_,_,_=fixture(audit['spec']);all_arrays={};quadrature=[]
    for order in orders:
        sample=TransportSampler(d,system.geometry.mesh,order);reference=mode.fields(sample.R,sample.Z)
        background=mode.ref.fields(sample.Z);p0=np.stack((background['electron_pressure'],background['ion_pressure']),axis=-1)
        p=p0*(d.thermal.gammas*reference['physical'][...,:1]/background['number_density'][...,None]+reference['physical'][...,4:6])
        targets=[*[reference['physical'][...,k] for k in range(9)],p[...,0],p[...,1],reference['current']]
        rows=[];fingerprints=[]
        for index in range(audit['steps']+1):
            states=[]
            for directory,record in zip(directories,records):
                filename='initial.npz' if index==0 else record['records'][index-1]['arrays_file']
                digest=record['initial_sha256'] if index==0 else record['records'][index-1]['arrays_sha256']
                if sha(directory/filename)!=digest:raise ValueError('Every original physical checkpoint required')
                states.append(dict(np.load(directory/filename)))
            minus,plus=states;when=index*audit['duration']/audit['steps'];growth=np.exp(mode.growth*when)
            fields,pressure,current=signal(sample,system,minus['state'],plus['state'],minus['magnetic'],plus['magnetic'],audit['epsilon'])
            values=[*[fields[...,k] for k in range(9)],pressure[...,0],pressure[...,1],current]
            row=[]
            for value,target in zip(values,targets):
                product=lambda a,b:np.sum(a*b,axis=-1) if a.ndim==3 else a*b
                denominator=float(np.sum(sample.weight*product(target,target)))
                gain=float(np.sum(sample.weight*product(value,target))/denominator)
                orthogonal=value-gain*target;error=value-growth*target
                row.append([gain,float(np.sqrt(np.sum(sample.weight*product(orthogonal,orthogonal))/denominator)),
                    float(np.sqrt(np.sum(sample.weight*product(error,error))/denominator))/growth])
            rows.append(row)
            fingerprints.append({name:fingerprint(value) for name,value in dict(fields=fields,pressure=pressure,current=current).items()})
        values=np.array(rows);all_arrays[f'q{order}_diagnostics']=values
        quadrature.append(dict(order=order,physical_fingerprints=fingerprints,
            maximum_relative_response_error={name:float(values[:,k,2].max()) for k,name in enumerate(NAMES)}))
        if order==16:
            endpoint=next(q for q in audit['quadrature'] if q['order']==16)['physical']['final']['errors']
            for k,name in enumerate(NAMES):
                if abs(values[-1,k,2]-endpoint[name]['relative_L2'])>1e-13:
                    raise RuntimeError('History endpoint disagrees with separately frozen physical audit')
    all_arrays['time']=np.arange(audit['steps']+1)*audit['duration']/audit['steps']
    all_arrays['expected_growth']=np.exp(mode.growth*all_arrays['time'])
    return dict(audit_directory=str(audit_directory),input_summary_sha256=sha(Path(audit_directory)/'summary.json'),
        spec=audit['spec'],steps=audit['steps'],epsilon=audit['epsilon'],duration=audit['duration'],quadrature=quadrature),all_arrays


def run(audits,destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    identity=provenance(['examples.matched_mode_history']);destination.mkdir(parents=True);freeze(destination,identity)
    cases=[]
    for audit in audits:
        record,arrays=calculate(audit);raw=destination/f'case_{len(cases):02}.npz';np.savez_compressed(raw,**arrays)
        record.update(arrays_file=raw.name,arrays_sha256=sha(raw));cases.append(record)
        print('physical response history',audit,flush=True)
    if provenance(['examples.matched_mode_history'])!=identity:raise RuntimeError('Source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],cases=cases,names=NAMES,
        columns=['projection_gain','orthogonal_error_over_initial_reference_norm','relative_error_to_growing_reference'])
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
