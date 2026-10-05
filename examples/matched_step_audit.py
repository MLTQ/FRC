"""Fresh original-equation substitution of saved performance-pilot stages."""
import json
from pathlib import Path
import numpy as np
from .matched_mode_fixture import fixture
from .matched_mode_rates import fingerprint
from .implicit_plasma_records import provenance, freeze, sha
from plasma_bench.modal_wall_coordinates import ModalWallCoordinates
from plasma_bench.uniform_wall_plasma import preconditioned


def run(directory,destination):
    directory=Path(directory);destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    record=json.loads((directory/'summary.json').read_text());source=json.loads((directory/'provenance.json').read_text())
    root=Path(__file__).resolve().parents[1]
    if not record['terminal']:raise ValueError('Terminal original pilot required')
    if any(sha(root/name)!=digest for name,digest in source['source_files'].items()):raise ValueError('Original source required')
    if sha(directory/'states.npz')!=record['arrays_sha256']:raise ValueError('Original arrays required')
    arrays=dict(np.load(directory/'states.npz'));mode,d,base,system,bias,_,_=fixture(record['spec'])
    system=preconditioned(system);c=ModalWallCoordinates(system,base,bias)
    state=arrays['initial_state'];magnetic=arrays['initial_magnetic'];system.derivative(state,c.field(magnetic))
    a=np.array([[.25,.25-np.sqrt(3)/6],[.25+np.sqrt(3)/6,.25]])
    identity=provenance(['examples.matched_step_audit']);destination.mkdir(parents=True);freeze(destination,identity)
    results=[]
    for index,original in enumerate(record['records']):
        actual=[];rates=[];residuals=[]
        for change in arrays[f'step_{index:03}_stages']:
            du,db=c.unpack(change);stage=state+du;field=c.field(magnetic+db)
            rate=system.derivative(stage,field)
            rates.append(c.pack(rate.fluid,rate.magnetic));residuals.append(rate.residual)
            actual.append(c.pack(stage-state,field.packed_increment(system.geometry)-magnetic))
        stage_error=np.array(actual)-record['dt']*a@rates
        next_state=arrays[f'step_{index:03}_state'];next_magnetic=arrays[f'step_{index:03}_magnetic']
        endpoint_error=c.pack(next_state-state,next_magnetic-magnetic)-.5*record['dt']*np.sum(rates,axis=0)
        maximum=max(float(np.max(abs(stage_error))),float(np.max(abs(endpoint_error))))
        if maximum>original['residual_limit']:raise RuntimeError('Saved actual stages fail original bound')
        results.append(dict(step=index+1,residual=maximum,limit=original['residual_limit'],
            induction_residuals=residuals,stage_error_fingerprint=fingerprint(stage_error),
            endpoint_error_fingerprint=fingerprint(endpoint_error)))
        state,magnetic=next_state,next_magnetic
    result=dict(terminal=True,input_directory=str(directory),input_summary_sha256=sha(directory/'summary.json'),
        input_source_sha256=source['source_sha256'],source_sha256=identity['source_sha256'],
        accepted_steps_verified=results,full_requested_pilot_complete=record['completed_steps']==record['requested_steps'])
    if provenance(['examples.matched_step_audit'])!=identity:raise RuntimeError('Audit source changed')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('pilot actual equations verified',directory,flush=True)
    return result
