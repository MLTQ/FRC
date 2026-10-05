"""Hold extracted physical fields fixed while subdividing the fluid mesh."""
import json
from pathlib import Path
import numpy as np
from plasma_bench.modal_wall_coordinates import ModalWallCoordinates
from plasma_bench.uniform_wall_plasma import preconditioned
from .matched_mode_fixture import fixture
from .radial_tangent_embedding import embed
from .current_transport_sampling import TransportSampler
from .extracted_response_fields import fields,inner
from .implicit_plasma_records import provenance,freeze,sha


def run(directory,destination):
    directory=Path(directory);destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    original=json.loads((directory/'summary.json').read_text());raw=directory/original['arrays_file']
    if not original['terminal'] or sha(raw)!=original['arrays_sha256']:raise ValueError('Complete original direction required')
    a=dict(np.load(raw));_,old,base0,system0,_,_,_=fixture(original['spec'])
    if not np.all(base0[...,0]==base0[0,0,:,0]):raise ValueError('Radially constant base density required')
    identity=provenance(['examples.extracted_fluid_refinement']);destination.mkdir(parents=True);freeze(destination,identity)
    cases=[]
    for elements in (4,8,16):
        spec=dict(original['spec'],elements=elements)
        _,d,base,system,bias,_,_=fixture(spec);system=preconditioned(system);c=ModalWallCoordinates(system,base,bias)
        du=embed(old,d,a['fluid_direction']);db=a['magnetic_direction'];rates=[];diagnostics=[]
        for epsilon in (1e-4,1e-5):
            pair=[system.derivative(base+sign*epsilon*du,c.field(sign*epsilon*db)) for sign in (-1,1)]
            rates.append(((pair[1].fluid-pair[0].fluid)/(2*epsilon),(pair[1].magnetic-pair[0].magnetic)/(2*epsilon)))
            diagnostics.append(dict(epsilon=epsilon,residuals=[p.residual for p in pair],induction=[p.diagnostics['induction_residuals'] for p in pair]))
        arrays=dict(fluid_direction=du,magnetic_direction=db,fluid_rate_eps1e4=rates[0][0],magnetic_rate_eps1e4=rates[0][1],
            fluid_rate_eps1e5=rates[1][0],magnetic_rate_eps1e5=rates[1][1]);quadrature=[]
        for order in (16,20):
            sample=TransportSampler(d,system.geometry.mesh,order);sample0=TransportSampler(old,system0.geometry.mesh,order)
            before=fields(sample0,system0,base0,a['fluid_direction'],db);value=fields(sample,system,base,du,db)
            tangents=[fields(sample,system,base,*rate) for rate in rates]
            baseline=fields(sample0,system0,base0,*ModalWallCoordinates(system0,base0,bias).unpack(a['rate_eps1e5']))
            channels={}
            for name,v in value.items():
                norm=inner(sample,v,v);diff=v-before[name]
                controls=[]
                for tangent in tangents:
                    rate=tangent[name];growth=inner(sample,v,rate)/norm;error=rate-growth*v
                    controls.append(dict(projected_growth=growth,off_direction_rate_relative=float(np.sqrt(inner(sample,error,error)/inner(sample,rate,rate)))))
                delta=tangents[1][name]-baseline[name]
                channels[name]=dict(embedding_relative_error=float(np.sqrt(inner(sample,diff,diff)/norm)),rhs_controls=controls,
                    rate_change_from_original=float(np.sqrt(inner(sample,delta,delta)/inner(sample,baseline[name],baseline[name]))))
            quadrature.append(dict(order=order,channels=channels))
        path=destination/f'case_{len(cases):02}.npz';np.savez_compressed(path,**arrays)
        cases.append(dict(spec=spec,quadrature=quadrature,rhs_diagnostics=diagnostics,arrays_file=path.name,arrays_sha256=sha(path)))
        print('complete same-field fluid refinement',elements,flush=True)
    if provenance(['examples.extracted_fluid_refinement'])!=identity:raise RuntimeError('Source changed')
    result=dict(terminal=True,source_sha256=identity['source_sha256'],input_summary_sha256=sha(directory/'summary.json'),cases=cases)
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
