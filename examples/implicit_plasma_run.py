"""Immutable individual implicit trajectory with accepted-state progress evidence."""
import argparse
import json
from pathlib import Path
import traceback
import numpy as np
from .implicit_plasma_evolution import run_case
from .implicit_plasma_records import provenance, freeze, sha


def run(directory, spec):
    directory = Path(directory)
    if directory.exists(): raise ValueError('Refusing to overwrite trajectory')
    seeds = [__name__ if __name__ != '__main__' else 'examples.implicit_plasma_run']
    identity = provenance(seeds); directory.mkdir(parents=True); freeze(directory, identity)
    latest = {}
    def checkpoint(history, arrays):
        latest.update(history)
        raw = directory/'progress.npz'; np.savez_compressed(raw, **arrays)
        (directory/'progress.json').write_text(json.dumps(dict(history, terminal=False,
            arrays_sha256=sha(raw), source_sha256=identity['source_sha256']), indent=2)+'\n')
        print('accepted', history['step'], 'residual', history['max_residual'], flush=True)
    try:
        result, arrays = run_case(spec, checkpoint)
        if provenance(seeds) != identity: raise RuntimeError('Source changed during trajectory')
        raw = directory/'case.npz'; np.savez_compressed(raw, **arrays)
        result.update(arrays_file=raw.name, arrays_sha256=sha(raw), source_sha256=identity['source_sha256'])
        (directory/'case.json').write_text(json.dumps(result, indent=2)+'\n')
    except Exception as exc:
        (directory/'failure.json').write_text(json.dumps(dict(spec=spec, status='failed', error=repr(exc),
            traceback=traceback.format_exc(), last_checkpoint=latest, source_sha256=identity['source_sha256']), indent=2)+'\n')
        raise
    print('complete', spec, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('directory')
    parser.add_argument('--kind', choices=('contact', 'wave'), required=True)
    parser.add_argument('--elements', type=int, required=True); parser.add_argument('--steps', type=int, default=120)
    parser.add_argument('--amplitude', type=float, default=1.)
    parser.add_argument('--method', choices=('midpoint', 'gauss'), default='midpoint')
    args = parser.parse_args(); run(args.directory, dict(kind=args.kind, elements=args.elements,
                                                       steps=args.steps, amplitude=args.amplitude, method=args.method))
