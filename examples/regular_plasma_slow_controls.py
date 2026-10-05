"""Frozen amplitude/mesh controls for the slowly evolving parent-wave branches."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from .regular_plasma_rates import run_case
from .regular_plasma_reference import initialize
from .regular_plasma_norms import measure
from .regular_plasma_records import provenance, freeze, sha


def run(destination):
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Refusing to overwrite slow-mode controls')
    identity = provenance(); identity['source_files']['examples/regular_plasma_slow_controls.py'] = sha(__file__)
    identity['source_sha256'] = hashlib.sha256(json.dumps(identity['source_files'], sort_keys=True).encode()).hexdigest()
    destination.mkdir(parents=True); freeze(destination, identity)
    records = []
    specs = [dict(kind='wave', elements=n, branch=b, amplitude=a)
             for b in (0, 1) for n, a in ((32, .01), (32, .001), (16, .001))]
    for i, spec in enumerate(specs):
        result, arrays = run_case(spec)
        path = destination/f'case_{i:02d}.npz'; np.savez_compressed(path, **arrays)
        repeated, replay = run_case(spec)
        assert all(np.array_equal(arrays[k], replay[k]) for k in arrays)
        assert all(result[k] == repeated[k] for k in result if k != 'elapsed_seconds')
        ref, d, u, system, field = initialize(spec); rate = system.derivative(u, field)
        result.update(source_sha256=identity['source_sha256'], arrays_file=path.name, arrays_sha256=sha(path),
                      bitwise_replayed=True, order16=measure(spec, ref, d, u, system, field, rate, 16),
                      order20=measure(spec, ref, d, u, system, field, rate, 20))
        path.with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n'); records.append(result)
        print(path.name, 'complete with replay and independent quadrature', flush=True)
    root = Path(__file__).resolve().parents[1]
    assert all(sha(root/name) == value for name, value in identity['source_files'].items())
    (destination/'summary.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'], terminal=True,
                                                          records=records), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('destination'); run(parser.parse_args().destination)
