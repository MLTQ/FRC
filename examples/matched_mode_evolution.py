"""Complete signed continuum-mode trajectories with accepted-step checkpoints."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from .matched_mode_fixture import fixture
from .implicit_plasma_records import provenance, freeze, sha
from plasma_bench.uniform_wall_plasma import preconditioned
from plasma_bench.recycled_modal_gauss import RecycledModalGauss


def run(destination, spec, *, epsilon=1e-4, steps=60, duration=1e-4, absolute_tolerance=1e-12):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    if not np.isfinite(epsilon) or epsilon == 0: raise ValueError('Finite nonzero signed amplitude required')
    if type(steps) is not int or steps < 1 or not np.isfinite(duration) or duration <= 0:
        raise ValueError('Positive steps and finite duration required')
    seeds = ['examples.matched_mode_evolution']; identity = provenance(seeds)
    destination.mkdir(parents=True); freeze(destination, identity)
    mode, _, base, system, bias, du, db = fixture(spec); system = preconditioned(system)
    solver = RecycledModalGauss(system, base, bias, absolute_tolerance=absolute_tolerance)
    state = base+epsilon*du; field = solver.coordinates.field(epsilon*db)
    initial = state.copy(); initial_b = field.packed_increment(system.geometry).copy()
    np.savez_compressed(destination/'initial.npz', state=state, magnetic=initial_b)
    record = dict(spec=spec, epsilon=epsilon, steps=steps, duration=duration, dt=duration/steps,
        absolute_tolerance=absolute_tolerance, relative_tolerance=solver.rtol,
        continuum_growth=mode.growth, input_sha256=mode.input_sha256,
        source_sha256=identity['source_sha256'], initial_sha256=sha(destination/'initial.npz'),
        completed_steps=0, records=[], terminal=False, success=False)
    (destination/'summary.json').write_text(json.dumps(record, indent=2)+'\n')
    for step in range(steps):
        start = time.monotonic()
        try:
            updated, evolved, diagnostic = solver.step(state, field, record['dt'])
        except Exception as error:
            record.update(failure_type=type(error).__name__, failure_message=str(error), failed_step=step+1)
            break
        stage_changes = diagnostic.pop('stage_changes')
        raw = destination/f'step_{step+1:04}.npz'
        np.savez_compressed(raw, state=updated, magnetic=evolved.packed_increment(system.geometry), stage_changes=stage_changes)
        state, field = updated, evolved
        record['completed_steps'] = step+1
        record['records'].append(dict(step=step+1, time=(step+1)*record['dt'],
            arrays_file=raw.name, arrays_sha256=sha(raw), elapsed_seconds=time.monotonic()-start, **diagnostic))
        (destination/'summary.json').write_text(json.dumps(record, indent=2)+'\n')
        if step == 0 or (step+1)%5 == 0:
            print(destination.name, step+1, '/', steps, 'RHS', diagnostic['rhs_evaluations'],
                'seconds', round(record['records'][-1]['elapsed_seconds'], 3), flush=True)
    if provenance(seeds) != identity: raise RuntimeError('Source changed during trajectory')
    np.savez_compressed(destination/'states.npz', initial_state=initial, initial_magnetic=initial_b,
        final_state=state, final_magnetic=field.packed_increment(system.geometry))
    record.update(terminal=True, success=record['completed_steps']==steps,
        final_time=record['completed_steps']*record['dt'], arrays_sha256=sha(destination/'states.npz'))
    (destination/'summary.json').write_text(json.dumps(record, indent=2)+'\n')
    print(destination.name, 'terminal', record['success'], record.get('failure_message', ''), flush=True)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination'); parser.add_argument('spec')
    parser.add_argument('--epsilon', type=float, default=1e-4); parser.add_argument('--steps', type=int, default=60)
    parser.add_argument('--duration', type=float, default=1e-4); parser.add_argument('--atol', type=float, default=1e-12)
    args = parser.parse_args()
    run(args.destination, json.loads(args.spec), epsilon=args.epsilon, steps=args.steps,
        duration=args.duration, absolute_tolerance=args.atol)
