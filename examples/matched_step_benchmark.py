"""Frozen full-equation nonlinear-step cost and failure records."""
import json
from pathlib import Path
import time
import numpy as np
from .matched_mode_fixture import fixture
from .implicit_plasma_records import provenance, freeze, sha
from plasma_bench.uniform_wall_plasma import preconditioned
from plasma_bench.modal_wall_gauss import ModalWallGauss
from plasma_bench.accelerated_modal_gauss import AcceleratedModalGauss


def run(destination, spec, *, epsilon=1e-4, intervals=60, accelerated=True, steps=1):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    identity = provenance(['examples.matched_step_benchmark'])
    destination.mkdir(parents=True); freeze(destination, identity)
    mode, d, base, system, bias, du, db = fixture(spec)
    system = preconditioned(system)
    solver_class = AcceleratedModalGauss if accelerated else ModalWallGauss
    solver = solver_class(system, base, bias); c = solver.coordinates
    state = base+epsilon*du; field = c.field(epsilon*db); dt = 1e-4/intervals
    records = []; initial = state.copy(); initial_b = field.packed_increment(system.geometry)
    arrays = dict(initial_state=initial, initial_magnetic=initial_b)
    result = dict(spec=spec, epsilon=epsilon, intervals=intervals, dt=dt,
        accelerated=accelerated, requested_steps=steps, continuum_growth=mode.growth,
        input_sha256=mode.input_sha256, source_sha256=identity['source_sha256'])
    for step in range(steps):
        start = time.monotonic()
        try:
            state, field, diagnostic = solver.step(state, field, dt)
        except Exception as error:
            result.update(failure_type=type(error).__name__, failure_message=str(error),
                failed_step=step, failed_step_seconds=time.monotonic()-start)
            break
        arrays[f'step_{step:03}_state'] = state
        arrays[f'step_{step:03}_magnetic'] = field.packed_increment(system.geometry)
        arrays[f'step_{step:03}_stages'] = diagnostic.pop('stage_changes')
        records.append(dict(step=step+1, elapsed_seconds=time.monotonic()-start, **diagnostic))
        print('step', step+1, 'seconds', records[-1]['elapsed_seconds'],
            'iterations', diagnostic['iterations'], 'RHS', diagnostic['rhs_evaluations'], flush=True)
    arrays.update(final_state=state, final_magnetic=field.packed_increment(system.geometry))
    np.savez_compressed(destination/'states.npz', **arrays)
    if provenance(['examples.matched_step_benchmark']) != identity: raise RuntimeError('Source changed')
    result.update(terminal=True, completed_steps=len(records), records=records, arrays_sha256=sha(destination/'states.npz'))
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print('terminal', len(records), '/', steps, result.get('failure_message', 'passed'), flush=True)
    return result
