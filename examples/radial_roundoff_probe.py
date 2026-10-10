"""Residual-path diagnosis on the exact rejected radial256 input."""
import ast
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.current_transport import CurrentTransport
from plasma_bench.regular_species import RegularSpecies
from plasma_bench.uniform_wall_plasma import UniformWallInduction
from plasma_bench.regular_plasma_step import field_from_packed
from plasma_bench.kernel_induction import separate
from .wall_roundoff_core import WallRoundoffCore
from .archived_contact_mode import ArchivedContactMode
from .radial_contact_fixture import fixture
from .implicit_plasma_records import provenance,freeze,sha


def run(destination,original='outputs/radial16_directional_refinement',index=1):
    destination=Path(destination);root=Path(original)
    if destination.exists():raise ValueError('Refusing overwrite')
    record=json.loads((root/f'case_{index:02}.json').read_text());path=root/record['arrays_file']
    if record['status']!='failed' or sha(path)!=record['arrays_sha256']:raise ValueError('Exact rejected input required')
    saved=dict(np.load(path));mode=ArchivedContactMode();identity=provenance(['examples.radial_roundoff_probe'])
    destination.mkdir(parents=True);freeze(destination,identity);start=time.monotonic()
    d,u,s,bias,du,db=fixture(record['spec'],mode)
    if any(not np.array_equal(saved[k],v) for k,v in [('state',u),('fluid_perturbation',du),('magnetic_perturbation',db)]):raise ValueError('Fixture mismatch')
    state=saved['minus_state'].copy();field=field_from_packed(s.geometry,saved['minus_magnetic'],bias.axial_bias)
    model=RegularCurrent(s.mapping,state[...,0],s.ion_mass);bulk,_=s.fluid.derivative(state)
    em=model.electromagnetic_pair(state,field,bulk[...,0]);tr=CurrentTransport(model).pair(state,em['current'],bulk[...,0]);sp=RegularSpecies(model).pair(state,em['current'])
    forcing=em['forcing']+tr['forcing']+sp['reduced_forcing'];solver=UniformWallInduction(model,**s.solver_options)
    electric0=solver.mass_inverse(forcing);small,kernel=separate(s.geometry,electric0);weak=model.mass@small
    result=dict(terminal=False,input_directory=str(root),input_arrays_sha256=sha(path),spec=record['spec'],source_sha256=identity['source_sha256'])
    try:
        accepted=WallRoundoffCore.solve(solver,weak)
        result.update(inner_status='accepted',inner_residual=accepted.residual)
    except RuntimeError as error:result.update(inner_status='rejected',inner_failure=str(error))
    attempts=ast.literal_eval(record['failure_message'].split('residual failure: ',1)[1].rsplit('; no input modified',1)[0])
    result['matches_original_first_attempt_failure']=result.get('inner_failure')==attempts[0].get('error')
    if not result['matches_original_first_attempt_failure']:raise ValueError('Instrumentation changed original failed attempt')
    arrays=dict(state=state,magnetic=saved['minus_magnetic'],forcing=forcing,weak_remainder=weak,**solver.roundoff_arrays)
    result.update(trace=solver.refinement_trace,roundoff=solver.roundoff_details,full_plasma_rate_accepted=False,
        unchanged_input=bool(np.array_equal(state,saved['minus_state'])),tolerance=s.solver_options['tolerance'])
    path=destination/'probe.npz';np.savez_compressed(path,**arrays)
    if provenance(['examples.radial_roundoff_probe'])!=identity:raise ValueError('Source changed')
    result.update(terminal=True,arrays_file=path.name,arrays_sha256=sha(path),elapsed_seconds=time.monotonic()-start)
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('roundoff probe',result['inner_status'],result['roundoff'],flush=True)
