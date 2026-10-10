"""Explicit quadrature, amplitude and residual controls for full signed rates."""
import json
from pathlib import Path
import numpy as np
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, cases, pairs=()):
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Refusing overwrite')
    hashes = {}; records = []; raw = []
    for directory, index in cases:
        directory = Path(directory); path = directory/'summary.json'
        summary = json.loads(path.read_text()); hashes[str(path)] = sha(path)
        if not summary['terminal']:
            raise ValueError('Terminal complete workflow required, including explicit failures')
        case = summary['cases'][index]
        source = json.loads((directory/'provenance.json').read_text())
        if summary['source_sha256'] != source['source_sha256'] or any(sha(p) != h for p, h in source['source_files'].items()):
            raise ValueError('Unchanged recorded source required')
        path = directory/case['arrays_file']
        if sha(path) != case['arrays_sha256']:
            raise ValueError('Original complete raw input required')
        hashes[str(path)] = sha(path); raw.append(dict(np.load(path)))
        record = dict(directory=str(directory), index=index, spec=case['spec'], status=case['status'],
            epsilon=case['epsilon'], reference_sha256=case['input_sha256'])
        if case['status'] == 'complete':
            bound = case['spec'].get('tolerance', 1e-10)
            if not np.isfinite(case['max_induction_residual']) or case['max_induction_residual'] > bound:
                raise ValueError('Original declared strict bound violated')
            record.update(max_induction_residual=case['max_induction_residual'], declared_bound=bound)
            low, high = case['quadrature']; differences = []
            for group in ('rate_errors', 'initial_errors'):
                for name in low[group]:
                    for region in ('relative_L2', 'axis_element', 'wall_element'):
                        a = low[group][name][region]; b = high[group][name][region]
                        differences.append(dict(group=group, name=name, region=region, low=a, high=b,
                            absolute=abs(a-b), relative=abs(a-b)/max(abs(a), abs(b), np.finfo(float).tiny)))
            record.update(quadrature_differences=differences,
                largest_relative=max(differences, key=lambda row: row['relative']),
                largest_absolute=max(differences, key=lambda row: row['absolute']),
                physical_rate_errors={k: v['relative_L2'] for k, v in high['rate_errors'].items()})
        elif case['status'] == 'failed':
            record.update(failed_sign=case['failed_sign'], failure_message=case['failure_message'])
        else:
            raise ValueError('Interrupted progress cannot pass a terminal control')
        records.append(record)
    comparisons = []; arrays = {}
    for label, first, second in pairs:
        a, b = records[first], records[second]
        keys = ('elements', 'points', 'field_nr', 'field_nz', 'quadrature')
        if a['reference_sha256'] != b['reference_sha256'] or any(a['spec'].get(k) != b['spec'].get(k) for k in keys):
            raise ValueError('Same physical grid/reference required for paired control')
        if a['status'] != 'complete' or b['status'] != 'complete':
            comparisons.append(dict(label=label, first=first, second=second, status='not_comparable_failed_solve'))
            continue
        differences = {}
        for name in ('fluid_rate', 'magnetic_rate', 'current_rate', 'electric'):
            x = (raw[first]['plus_'+name]-raw[first]['minus_'+name])/(2*a['epsilon'])
            y = (raw[second]['plus_'+name]-raw[second]['minus_'+name])/(2*b['epsilon'])
            arrays[label+'_'+name+'_difference'] = x-y
            differences[name] = float(np.linalg.norm(x-y)/max(np.linalg.norm(y), np.finfo(float).tiny))
        comparisons.append(dict(label=label, first=first, second=second, status='complete',
            coefficient_relative_differences=differences,
            physical_error_norm_differences={name: b['physical_rate_errors'][name]-a['physical_rate_errors'][name]
                for name in a['physical_rate_errors']}))
    identity = provenance(['examples.joint_rate_controls'])
    destination.mkdir(parents=True); freeze(destination, identity)
    path = destination/'differences.npz'; np.savez_compressed(path, **arrays)
    result = dict(terminal=True, source_sha256=identity['source_sha256'], input_hashes=hashes,
        cases=records, comparisons=comparisons, arrays_file=path.name, arrays_sha256=sha(path),
        coefficient_metric='Euclidean coefficient norms; distinct from independent physical field error norms')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
