"""Complete exact-input observation of the current split-inverse rejection."""
import ast
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.current_transport import CurrentTransport
from plasma_bench.regular_species import RegularSpecies
from plasma_bench.direct_metric_wall_plasma import DirectMetricWallInduction
from plasma_bench.regular_plasma_step import field_from_packed
from plasma_bench.kernel_induction import separate
from .split_residual_observer import SplitResidualObserver
from .archived_contact_mode import ArchivedContactMode
from .radial_contact_fixture import fixture
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, original='outputs/radial16_finest_controls', index=1):
    destination = Path(destination); root = Path(original)
    if destination.exists(): raise ValueError('Refusing overwrite')
    summary = json.loads((root/'summary.json').read_text())
    record = summary['cases'][index]
    source = json.loads((root/'provenance.json').read_text())
    if not summary['terminal'] or record['status'] != 'failed' or record['failed_sign'] != 'minus':
        raise ValueError('Original terminal minus rejection required')
    if summary['source_sha256'] != source['source_sha256'] or any(sha(p) != h for p, h in source['source_files'].items()):
        raise ValueError('Original source closure required')
    path = root/record['arrays_file']
    if sha(path) != record['arrays_sha256']: raise ValueError('Exact rejected input required')
    saved = dict(np.load(path))
    mode = ArchivedContactMode(summary['reference_directory'], summary['reference_index'])
    if mode.input_sha256 != record['input_sha256'] or mode.input_summary_sha256 != record['input_summary_sha256']:
        raise ValueError('Same independent continuum reference required')
    identity = provenance(['examples.split_residual_probe'])
    destination.mkdir(parents=True); freeze(destination, identity); start = time.monotonic()
    result = dict(terminal=False, input_directory=str(root), input_arrays_sha256=sha(path),
        input_summary_sha256=sha(root/'summary.json'), spec=record['spec'], source_sha256=identity['source_sha256'])
    def save(): (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    save(); print('Exact split residual probe initializing', flush=True)
    d, u, system, bias, du, db = fixture(record['spec'], mode)
    exact = dict(state=u, fluid_perturbation=du, magnetic_perturbation=db,
        minus_state=u-record['epsilon']*du, minus_magnetic=-record['epsilon']*db)
    if any(not np.array_equal(saved[k], v) for k, v in exact.items()):
        raise ValueError('Bitwise fixture/signed-input identity required')
    state = saved['minus_state'].copy(); magnetic = saved['minus_magnetic'].copy()
    field = field_from_packed(system.geometry, magnetic, bias.axial_bias)
    model = RegularCurrent(system.mapping, state[..., 0], system.ion_mass)
    bulk, _ = system.fluid.derivative(state)
    em = model.electromagnetic_pair(state, field, bulk[..., 0])
    tr = CurrentTransport(model).pair(state, em['current'], bulk[..., 0])
    sp = RegularSpecies(model).pair(state, em['current'])
    forcing = em['forcing']+tr['forcing']+sp['reduced_forcing']
    solver = DirectMetricWallInduction(model, **system.solver_options)
    electric0 = solver.mass_inverse(forcing); small, kernel = separate(system.geometry, electric0)
    weak = model.mass@small
    print('Exact split residual probe observing original first attempt', flush=True)
    try:
        accepted = SplitResidualObserver.solve(solver, weak)
        result.update(inner_status='accepted', inner_residual=accepted.residual)
    except RuntimeError as error:
        result.update(inner_status='rejected', inner_failure=str(error))
    attempts = ast.literal_eval(record['failure_message'].split('residual failure: ', 1)[1].rsplit('; no input modified', 1)[0])
    matches = result.get('inner_failure') == attempts[0].get('error')
    unchanged = np.array_equal(state, saved['minus_state']) and np.array_equal(magnetic, saved['minus_magnetic'])
    arrays = dict(state=state, magnetic=magnetic, forcing=forcing, weak_remainder=weak,
        separated_kernel=kernel, **solver.observation_arrays)
    path = destination/'probe.npz'; np.savez_compressed(path, **arrays)
    if provenance(['examples.split_residual_probe']) != identity: raise ValueError('Source changed')
    result.update(terminal=True, observation_consistent=bool(matches and unchanged),
        matches_original_first_attempt_failure=matches, unchanged_input=bool(unchanged),
        full_plasma_rate_accepted=False, trace=solver.refinement_trace, observations=solver.observations,
        arrays_file=path.name, arrays_sha256=sha(path), elapsed_seconds=time.monotonic()-start)
    save()
    if not result['observation_consistent']: raise ValueError('Observation did not preserve original rejected workflow')
    print('Exact split residual probe terminal; original rejection matched', flush=True)
    return result
