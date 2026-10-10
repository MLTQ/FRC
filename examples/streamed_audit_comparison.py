"""Measure diagnostic batching differences with identical full physical rates."""
import json
from pathlib import Path
import numpy as np
from .implicit_plasma_records import provenance, freeze, sha


def load(directory, index):
    directory = Path(directory)
    summary = json.loads((directory/'summary.json').read_text())
    if not summary['terminal']:
        raise ValueError('Complete full workflow required')
    case = summary['cases'][index]
    if case['status'] != 'complete':
        raise ValueError('Successful full signed case required')
    source = json.loads((directory/'provenance.json').read_text())
    if summary['source_sha256'] != source['source_sha256']:
        raise ValueError('Matching source closure required')
    if any(sha(name) != digest for name, digest in source['source_files'].items()):
        raise ValueError('Unchanged original source required')
    path = directory/case['arrays_file']
    if sha(path) != case['arrays_sha256']:
        raise ValueError('Unchanged full arrays required')
    return case, dict(np.load(path))


def run(original, repeated, destination, original_index=0, repeated_index=0):
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Refusing overwrite')
    before, a = load(original, original_index); after, b = load(repeated, repeated_index)
    for key in ('spec', 'epsilon', 'input_sha256', 'input_summary_sha256', 'radial_mode',
                'reference_points', 'continuum_growth', 'rate_diagnostics', 'max_induction_residual'):
        if before[key] != after[key]:
            raise ValueError('Identical physical solve required: '+key)
    if a.keys() != b.keys() or any(not np.array_equal(value, b[key]) for key, value in a.items()):
        raise ValueError('All initial/signed states, full rates and face-integral arrays must agree bitwise')
    metrics = {}; fingerprints = {}
    def compare(first, second, path):
        if isinstance(first, dict):
            if first.keys() != second.keys():
                raise ValueError('Identical diagnostic schema required: '+path)
            for key in first:
                compare(first[key], second[key], path+'/'+key)
        elif isinstance(first, list):
            if len(first) != len(second):
                raise ValueError('Identical quadrature points required')
            for index, (left, right) in enumerate(zip(first, second)):
                compare(left, right, path+'/'+str(index))
        elif isinstance(first, str):
            fingerprints[path] = first == second
        elif isinstance(first, (int, float)):
            if not np.isfinite([first, second]).all():
                raise ValueError('Finite physical comparison required')
            absolute = abs(first-second)
            metrics[path] = dict(original=first, streamed=second, absolute_difference=absolute,
                relative_difference=absolute/max(abs(first), abs(second), np.finfo(float).tiny))
        elif first != second:
            raise ValueError('Unexpected diagnostic type: '+path)
    compare(before['quadrature'], after['quadrature'], 'quadrature')
    identity = provenance(['examples.streamed_audit_comparison'])
    destination.mkdir(parents=True); freeze(destination, identity)
    result = dict(terminal=True, source_sha256=identity['source_sha256'],
        original_directory=str(original), repeated_directory=str(repeated),
        original_summary_sha256=sha(Path(original)/'summary.json'),
        repeated_summary_sha256=sha(Path(repeated)/'summary.json'),
        arrays_bitwise_equal=True, solver_diagnostics_bitwise_equal=True,
        numeric_differences=metrics, fingerprints_equal=fingerprints,
        maximum_absolute_difference=max(v['absolute_difference'] for v in metrics.values()),
        maximum_relative_difference=max(v['relative_difference'] for v in metrics.values()))
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
