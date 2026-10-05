"""Independent axis/wall cell-budget norms from complete physical flux archives."""
import json
from pathlib import Path
import numpy as np
from .matched_mode_fixture import fixture
from .current_transport_sampling import TransportSampler
from .physical_cell_budget import NAMES
from .implicit_plasma_records import provenance,freeze,sha


def run(directory,destination):
    directory=Path(directory);destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    original=json.loads((directory/'summary.json').read_text())
    if not original['terminal']:raise ValueError('Complete interface study required')
    identity=provenance(['examples.matched_interface_regions']);destination.mkdir(parents=True);freeze(destination,identity)
    cases=[]
    for record in original['cases']:
        path=directory/record['arrays_file']
        if sha(path)!=record['arrays_sha256']:raise ValueError('Original cell arrays required')
        arrays=dict(np.load(path));_,d,_,system,_,_,_=fixture(record['spec']);nr=system.geometry.mesh.nr
        width=nr//d.elements;quadrature=[]
        for parent in record['quadrature']:
            order=parent['order'];sample=TransportSampler(d,system.geometry.mesh,order)
            _,gradient=sample.electric(arrays['current'],True)
            gradient_cell=(sample.weight*np.sum(gradient*gradient,axis=(-1,-2))).reshape(nr,order,-1,order).sum(axis=(1,3))
            volume=arrays[f'q{order}_volume'];residual=arrays[f'q{order}_cell_balance']
            rate=arrays[f'q{order}_volume_rate'];out=arrays[f'q{order}_cell_outflux'];regions={}
            for name,mask in [('axis',np.arange(nr)<width),('wall',np.arange(nr)>=nr-width)]:
                values={}
                for k,channel in enumerate(NAMES):
                    numerator=float(np.sum(residual[mask,:,k]**2/volume[mask]))
                    denominator=float(gradient_cell[mask].sum()) if channel=='current' else float(np.sum((rate[mask,:,k]**2+out[mask,:,k]**2)/volume[mask]))
                    values[channel]=dict(relative_L2=float(np.sqrt(numerator/max(denominator,np.finfo(float).tiny))),
                        absolute_L2=float(np.sqrt(numerator)),maximum_density=float(np.max(abs(residual[mask,:,k]/volume[mask]))))
                regions[name]=values
            quadrature.append(dict(order=order,regions=regions))
        cases.append(dict(spec=record['spec'],epsilon=record['epsilon'],arrays_sha256=sha(path),quadrature=quadrature))
    result=dict(terminal=True,source_sha256=identity['source_sha256'],input_source_sha256=original['source_sha256'],
        input_summary_sha256=sha(directory/'summary.json'),cases=cases)
    if provenance(['examples.matched_interface_regions'])!=identity:raise RuntimeError('Source changed')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
