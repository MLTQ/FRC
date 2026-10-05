"""Physical time/amplitude differences between audited signed trajectories."""
import json
from pathlib import Path
import numpy as np
from .matched_mode_fixture import fixture
from .matched_mode_physical import signal
from .matched_mode_rates import fingerprint
from .current_transport_sampling import TransportSampler
from .regular_plasma_norms import norm
from .implicit_plasma_records import provenance, freeze, sha
from plasma_bench.regular_plasma_step import field_from_packed


def load(audit_directory):
    audit_directory=Path(audit_directory);record=json.loads((audit_directory/'summary.json').read_text())
    if not record['terminal']:raise ValueError('Complete independent audit required')
    arrays=[]
    for directory, digest in zip(record['input_directories'],record['input_summary_sha256']):
        directory=Path(directory)
        if sha(directory/'summary.json')!=digest:raise ValueError('Unchanged audited trajectory required')
        summary=json.loads((directory/'summary.json').read_text())
        if sha(directory/'states.npz')!=summary['arrays_sha256']:raise ValueError('Unchanged raw states required')
        arrays.append(dict(np.load(directory/'states.npz')))
    return record,arrays


def physical(sample,system,bias,arrays,epsilon):
    minus,plus=arrays
    fields,pressure,current=signal(sample,system,minus['final_state'],plus['final_state'],
        minus['final_magnetic'],plus['final_magnetic'],epsilon)
    rates=[system.derivative(a['final_state'],field_from_packed(system.geometry,a['final_magnetic'],bias.axial_bias)) for a in arrays]
    electric=sample.electric((rates[1].electric-rates[0].electric)/(2*epsilon))
    return dict(fields=fields,pressure=pressure,current=current,electric=electric)


def run(first,second,destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    a,aa=load(first);b,ba=load(second)
    if a['spec']!=b['spec'] or a['duration']!=b['duration']:raise ValueError('Same meshes and duration required')
    identity=provenance(['examples.matched_mode_comparison']);destination.mkdir(parents=True);freeze(destination,identity)
    mode,d,_,system,bias,_,_=fixture(a['spec']);results=[]
    names=['density','radial_velocity','azimuthal_velocity','axial_velocity','electron_entropy','ion_entropy','radial_B','azimuthal_B','axial_B']
    for order in (16,20):
        sample=TransportSampler(d,system.geometry.mesh,order)
        values=[physical(sample,system,bias,arrays,record['epsilon']) for record,arrays in ((a,aa),(b,ba))]
        differences={name:norm(sample,values[1]['fields'][...,k],values[0]['fields'][...,k]) for k,name in enumerate(names)}
        differences.update({name:norm(sample,values[1][name],values[0][name]) for name in ('current','electric')})
        for k,name in enumerate(('electron_pressure','ion_pressure')):
            differences[name]=norm(sample,values[1]['pressure'][...,k],values[0]['pressure'][...,k])
        results.append(dict(order=order,differences=differences,
            physical_fingerprints=[{name:fingerprint(value) for name,value in values.items()} for values in values]))
    result=dict(terminal=True,source_sha256=identity['source_sha256'],first=str(first),second=str(second),
        input_sha256=[sha(Path(p)/'summary.json') for p in (first,second)],spec=a['spec'],duration=a['duration'],
        steps=[a['steps'],b['steps']],epsilon=[a['epsilon'],b['epsilon']],quadrature=results)
    if provenance(['examples.matched_mode_comparison'])!=identity:raise RuntimeError('Source changed')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('complete physical comparison',destination,flush=True)
    return result
