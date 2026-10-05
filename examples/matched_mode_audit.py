"""Independent signed trajectory comparison with original continuum growth."""
import argparse
import json
from pathlib import Path
import numpy as np
from .matched_mode_fixture import fixture
from .matched_mode_physical import changes, signal
from .matched_mode_rates import fingerprint
from .current_transport_sampling import TransportSampler
from .regular_plasma_norms import norm
from .implicit_plasma_records import provenance, freeze, sha
from plasma_bench.regular_plasma_step import field_from_packed


def measure(spec, negative, positive, epsilon, time, order):
    mode, d, _, system, bias, _, _ = fixture(spec)
    sample = TransportSampler(d, system.geometry.mesh, order)
    target = mode.fields(sample.R, sample.Z); factor = np.exp(mode.growth*time)
    names = ['density','radial_velocity','azimuthal_velocity','axial_velocity',
        'electron_entropy','ion_entropy','radial_B','azimuthal_B','axial_B']
    records = {}; physical_hashes = {}
    for prefix, growth in [('initial',1.),('final',factor)]:
        uminus = negative[prefix+'_state']; uplus = positive[prefix+'_state']
        bminus = negative[prefix+'_magnetic']; bplus = positive[prefix+'_magnetic']
        physical, pressure, current = signal(sample, system, uminus, uplus, bminus, bplus, epsilon)
        expected = growth*target['physical']
        errors = {name:norm(sample, physical[..., k], expected[..., k]) for k,name in enumerate(names)}
        errors['current'] = norm(sample,current,growth*target['current'])
        rates=[]
        for u,b in [(uminus,bminus),(uplus,bplus)]:
            rates.append(system.derivative(u,field_from_packed(system.geometry,b,bias.axial_bias)))
        electric=sample.electric((rates[1].electric-rates[0].electric)/(2*epsilon))
        errors['electric']=norm(sample,electric,growth*target['electric'])
        background=mode.ref.fields(sample.Z)
        p0=np.stack((background['electron_pressure'],background['ion_pressure']),axis=-1)
        expected_pressure=growth*p0*(d.thermal.gammas*target['physical'][..., :1]/background['number_density'][...,None]+target['physical'][...,4:6])
        for k,name in enumerate(('electron_pressure','ion_pressure')):
            errors[name]=norm(sample,pressure[...,k],expected_pressure[...,k])
        gains={name:float(np.sum(sample.weight*physical[...,k]*target['physical'][...,k])/
            np.sum(sample.weight*target['physical'][...,k]**2)) for k,name in enumerate(names)}
        records[prefix]=dict(errors=errors, projected_growth_factors=gains,
            expected_growth=growth, induction_residuals=[r.residual for r in rates],
            full_rate_diagnostics=[r.diagnostics for r in rates])
        physical_hashes[prefix]={name:fingerprint(value) for name,value in
            dict(physical=physical,pressure=pressure,current=current,electric=electric,expected=expected).items()}
    conservation={}
    for name,arrays in [('minus',negative),('plus',positive)]:
        conservation[name]=changes(sample,system,bias.axial_bias,arrays['initial_state'],arrays['initial_magnetic'],
            arrays['final_state'],arrays['final_magnetic'])
    return dict(order=order,physical=records,conservation=conservation,field_fingerprints=physical_hashes)


def run(negative_directory,positive_directory,destination):
    negative_directory=Path(negative_directory);positive_directory=Path(positive_directory);destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    records=[]; arrays=[]; sources=[];root=Path(__file__).resolve().parents[1]
    for directory in (negative_directory,positive_directory):
        record=json.loads((directory/'summary.json').read_text());source=json.loads((directory/'provenance.json').read_text())
        if not record['terminal'] or not record['success']:raise ValueError('Complete successful trajectories required')
        if any(sha(root/name)!=digest for name,digest in source['source_files'].items()):raise ValueError('Original source required')
        if sha(directory/'states.npz')!=record['arrays_sha256']:raise ValueError('Original raw arrays required')
        for step in record['records']:
            if sha(directory/step['arrays_file'])!=step['arrays_sha256']:raise ValueError('Complete step archives required')
        records.append(record);arrays.append(dict(np.load(directory/'states.npz')));sources.append(source)
    minus,plus=records
    for key in ('spec','steps','duration','absolute_tolerance','relative_tolerance','input_sha256','source_sha256'):
        if minus[key]!=plus[key]:raise ValueError('Matching signed inputs required: '+key)
    if minus['epsilon']!=-plus['epsilon'] or plus['epsilon']<=0:raise ValueError('Opposite equal amplitudes required')
    seeds=['examples.matched_mode_audit'];identity=provenance(seeds)
    if any(identity[k]!=sources[0][k] for k in ('python','numpy','scipy')):raise ValueError('Original runtime required')
    destination.mkdir(parents=True);freeze(destination,identity)
    result=dict(input_directories=[str(negative_directory),str(positive_directory)],
        input_source_sha256=minus['source_sha256'],source_sha256=identity['source_sha256'],
        input_summary_sha256=[sha(directory/'summary.json') for directory in (negative_directory,positive_directory)],
        spec=minus['spec'],epsilon=plus['epsilon'],duration=plus['duration'],steps=plus['steps'],
        quadrature=[measure(plus['spec'],*arrays,plus['epsilon'],plus['duration'],order) for order in (16,20)])
    if provenance(seeds)!=identity:raise RuntimeError('Audit source changed')
    result['terminal']=True
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('matched independent audit complete',destination,flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('negative');parser.add_argument('positive');parser.add_argument('destination')
    args=parser.parse_args();run(args.negative,args.positive,args.destination)
