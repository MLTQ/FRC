"""Archive bitwise comparisons of independently repeated wall diagnostics."""
import json
from pathlib import Path
import numpy as np
from .contact_wall_diagnosis import run as wall_run
from .continuum_wall_probe import run as continuum_run
from .implicit_plasma_audit import without_runtime
from .implicit_plasma_records import provenance, freeze, sha


def verify(original, repeated, destination):
    original = Path(original); repeated = Path(repeated); destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    before = json.loads((original/'summary.json').read_text())
    after = json.loads((repeated/'summary.json').read_text())
    source = json.loads((original/'provenance.json').read_text())
    if source != json.loads((repeated/'provenance.json').read_text()): raise ValueError('Same source/runtime required')
    root = Path(__file__).resolve().parents[1]
    if any(sha(root/name) != digest for name, digest in source['source_files'].items()): raise ValueError('Original source closure required')
    if not before['terminal'] or not after['terminal'] or without_runtime(before) != without_runtime(after):
        raise ValueError('Complete bitwise non-time metrics required')
    names = sorted(p.name for p in original.glob('*.npz'))
    if names != sorted(p.name for p in repeated.glob('*.npz')): raise ValueError('Complete matching raw files required')
    arrays = {}
    for name in names:
        a = dict(np.load(original/name)); b = dict(np.load(repeated/name))
        if a.keys() != b.keys() or any(not np.array_equal(value, b[key]) for key, value in a.items()):
            raise RuntimeError(f'Array replay failed: {name}')
        arrays[name] = dict(bitwise_equal=True, original_sha256=sha(original/name), repeated_sha256=sha(repeated/name))
    identity = provenance(['examples.wall_diagnostic_replay']); destination.mkdir(parents=True); freeze(destination, identity)
    result = dict(terminal=True, bitwise_metrics=True, arrays=arrays, original_directory=str(original),
        repeated_directory=str(repeated), input_source_sha256=source['source_sha256'],
        original_summary_sha256=sha(original/'summary.json'), repeated_summary_sha256=sha(repeated/'summary.json'),
        source_sha256=identity['source_sha256'])
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')


def replay(kind, original, repeated, destination):
    functions = {'wall': wall_run, 'continuum': continuum_run}
    if kind not in functions: raise ValueError('Known wall diagnostic required')
    functions[kind](repeated)
    verify(original, repeated, destination)


if __name__ == '__main__':
    import sys
    verify(*sys.argv[1:])
