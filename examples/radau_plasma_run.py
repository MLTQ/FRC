"""Immutable Radau trajectories, complete replay and independent quadrature."""
import argparse
import json
from pathlib import Path
import traceback
import numpy as np
from .radau_plasma_evolution import run_case
from .implicit_plasma_records import provenance, freeze, sha
from .implicit_plasma_audit import measure, without_runtime


SEEDS = ['examples.radau_plasma_run']


def run(directory, spec):
    directory = Path(directory)
    if directory.exists(): raise ValueError('Refusing to overwrite trajectory')
    identity = provenance(SEEDS); directory.mkdir(parents=True); freeze(directory, identity)
    latest = {}
    def checkpoint(history, arrays):
        latest.update(history); raw = directory/'progress.npz'; np.savez_compressed(raw, **arrays)
        (directory/'progress.json').write_text(json.dumps(dict(history, terminal=False,
            arrays_sha256=sha(raw), source_sha256=identity['source_sha256']), indent=2)+'\n')
        print('accepted', history['step'], 'residual', history['max_residual'], flush=True)
    try:
        result, arrays = run_case(spec, checkpoint)
        if provenance(SEEDS) != identity: raise RuntimeError('Source changed during trajectory')
        raw = directory/'case.npz'; np.savez_compressed(raw, **arrays)
        result.update(arrays_file=raw.name, arrays_sha256=sha(raw), source_sha256=identity['source_sha256'])
        (directory/'case.json').write_text(json.dumps(result, indent=2)+'\n')
    except Exception as exc:
        (directory/'failure.json').write_text(json.dumps(dict(spec=spec, status='failed', error=repr(exc),
            traceback=traceback.format_exc(), last_checkpoint=latest, source_sha256=identity['source_sha256']), indent=2)+'\n')
        raise
    print('complete', spec, flush=True)


def audit(directory, destination, replay=False):
    directory = Path(directory); destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing to overwrite audit')
    identity = provenance(SEEDS); original = json.loads((directory/'provenance.json').read_text())
    if identity != original: raise ValueError('Original source/runtime required')
    record = json.loads((directory/'case.json').read_text()); raw = directory/record['arrays_file']
    if sha(raw) != record['arrays_sha256']: raise ValueError('Original arrays hash required')
    arrays = dict(np.load(raw)); destination.mkdir(parents=True); freeze(destination, identity)
    result = dict(input_directory=str(directory), source_sha256=identity['source_sha256'],
        arrays_sha256=sha(raw), spec=record['spec'],
        order16=measure(record['spec'], arrays, record['duration'], 16),
        order20=measure(record['spec'], arrays, record['duration'], 20), full_replay=False)
    if replay:
        repeated, repeated_arrays = run_case(record['spec'])
        same = {k: np.array_equal(v, repeated_arrays[k]) for k, v in arrays.items()}
        metrics = without_runtime(repeated) == without_runtime({k: record[k] for k in repeated})
        result.update(full_replay=True, bitwise_arrays=same, bitwise_metrics=metrics)
        if not all(same.values()) or not metrics: raise RuntimeError('Complete trajectory replay differed')
    if provenance(SEEDS) != identity: raise RuntimeError('Audit source changed')
    result['terminal'] = True
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print('independent Radau audit complete', directory, 'replay', replay, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); actions = parser.add_subparsers(dest='action', required=True)
    p = actions.add_parser('run'); p.add_argument('directory')
    p.add_argument('--kind', choices=('contact', 'wave'), required=True)
    p.add_argument('--elements', type=int, required=True); p.add_argument('--steps', type=int, default=120)
    p.add_argument('--amplitude', type=float, default=1.)
    p = actions.add_parser('audit'); p.add_argument('directory'); p.add_argument('destination')
    p.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    if args.action == 'run':
        run(args.directory, dict(kind=args.kind, elements=args.elements, steps=args.steps, amplitude=args.amplitude, method='radau'))
    else: audit(args.directory, args.destination, args.replay)
