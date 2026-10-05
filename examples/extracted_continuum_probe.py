"""Independent continuum action on the physical extracted response and rate."""
import json
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0
from .matched_mode_fixture import fixture
from .current_transport_sampling import TransportSampler
from .extracted_response_fields import fields, inner
from .contact_mode_fields import magnetic_curl
from .contact_mode_projection import project
from .implicit_plasma_records import provenance, freeze, sha
from plasma_bench.modal_wall_coordinates import ModalWallCoordinates


def run(directory, destination):
    directory=Path(directory);destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    record=json.loads((directory/'summary.json').read_text());raw=directory/record['arrays_file']
    if not record['terminal'] or sha(raw)!=record['arrays_sha256']:raise ValueError('Complete original direction required')
    a=dict(np.load(raw));_,d,base,system,bias,_,_=fixture(record['spec']);c=ModalWallCoordinates(system,base,bias)
    identity=provenance(['examples.extracted_continuum_probe']);destination.mkdir(parents=True);freeze(destination,identity)
    cases=[]
    for order,points,modes in ((16,32,8),(16,64,16),(20,64,16),(20,64,32)):
        sample=TransportSampler(d,system.geometry.mesh,order)
        direction=fields(sample,system,base,a['fluid_direction'],a['magnetic_direction'])
        rate=fields(sample,system,base,*c.unpack(a['rate_eps1e5']))
        physical=lambda f:np.concatenate((f['density'],f['velocity'],f['entropy'],f['magnetic']),axis=-1)
        common=dict(radius=sample.r,axial=sample.z,radial_weight=sample.weight.sum(axis=1)/d.length,
            axial_weight=sample.weight.sum(axis=0)/(np.pi*d.radius**2))
        p=project(dict(common,physical_mode=physical(direction)),points,modes)
        r=project(dict(common,physical_mode=physical(rate)),points,modes)
        weights=p['radial_norm'][:,None,None]*d.length/points
        integral=lambda x,y:np.sum(weights*x*y,axis=(0,2))
        norm=integral(p['projected'],p['projected']);rate_norm=integral(r['projected'],r['projected'])
        predicted_norm=integral(p['predicted'],p['predicted'])
        difference=r['projected']-p['predicted']
        observed_growth=integral(p['projected'],r['projected'])/norm
        continuum_growth=integral(p['projected'],p['predicted'])/norm
        weak=direction['current'];strong=magnetic_curl(sample,a['magnetic_direction'])/MU0
        result=dict(order=order,points=points,radial_modes=modes,
            captured_direction_fractions=(p['captured_norms']/p['original_norms']).tolist(),
            captured_rate_fractions=(r['captured_norms']/r['original_norms']).tolist(),
            projected_actual_growth=observed_growth.tolist(),projected_continuum_growth=continuum_growth.tolist(),
            actual_continuum_rate_relative_difference=np.sqrt(integral(difference,difference)/rate_norm).tolist(),
            continuum_rate_norm_over_actual=np.sqrt(predicted_norm/rate_norm).tolist(),
            radial_mode_power_fractions=(np.sum(weights*p['projected']**2,axis=2)/norm).tolist(),
            strong_ampere_relative=float(np.sqrt(inner(sample,weak-strong,weak-strong)/inner(sample,weak,weak))))
        arrays={f'direction_{k}':v for k,v in p.items()};arrays.update({f'rate_{k}':v for k,v in r.items()})
        path=destination/f'case_{len(cases):02}.npz';np.savez_compressed(path,**arrays)
        result.update(arrays_file=path.name,arrays_sha256=sha(path));cases.append(result)
        print('complete continuum projection',order,points,modes,flush=True)
    if provenance(['examples.extracted_continuum_probe'])!=identity:raise RuntimeError('Source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],input_directory=str(directory),
        input_summary_sha256=sha(directory/'summary.json'),input_arrays_sha256=sha(raw),cases=cases,
        component_order=['density','radial_velocity','azimuthal_velocity','axial_velocity','electron_entropy','ion_entropy','radial_B','azimuthal_B','axial_B'])
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
