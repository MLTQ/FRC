"""Baseline, amplitude and quadrature controls for split-inverse full rates."""
import json
from pathlib import Path
import numpy as np
from .implicit_plasma_records import provenance, freeze, sha


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    inputs = {}
    def load(directory, index):
        root = Path('outputs')/directory; path = root/'summary.json'
        summary = json.loads(path.read_text()); inputs[str(path)] = sha(path)
        if not summary['terminal']: raise ValueError('Terminal batch required')
        case = summary['cases'][index]; path = root/case['arrays_file']
        if case['status'] != 'complete' or sha(path) != case['arrays_sha256']:
            raise ValueError('Hash-verified accepted full rates required')
        inputs[str(path)] = sha(path)
        return case, dict(np.load(path))
    old, oldraw = load('radial16_balanced_rates', 1)
    base, baseraw = load('radial16_split_inverse_rates', 0)
    fine, fineraw = load('radial16_split_inverse_rates', 1)
    small, smallraw = load('radial16_split_inverse_rates', 2)
    fixed = ['state', 'fluid_perturbation', 'magnetic_perturbation', 'minus_state', 'plus_state',
             'minus_magnetic', 'plus_magnetic', 'minus_fluid_rate', 'plus_fluid_rate']
    if old['spec'] != base['spec'] or any(not np.array_equal(oldraw[k], baseraw[k]) for k in fixed):
        raise ValueError('Actual baseline inputs or original fluid rate changed')
    if any(fine['spec'][k] != small['spec'][k] for k in ('elements', 'points', 'field_nr', 'field_nz')):
        raise ValueError('Same mesh amplitude control required')
    arrays = {}; controls = {}
    for label, first, a, second, b in [('baseline', old, oldraw, base, baseraw),
                                      ('amplitude', fine, fineraw, small, smallraw)]:
        rows = {}
        for key in ('fluid_rate', 'magnetic_rate', 'current_rate', 'electric'):
            x = (a['plus_'+key]-a['minus_'+key])/(2*first['epsilon'])
            y = (b['plus_'+key]-b['minus_'+key])/(2*second['epsilon'])
            arrays[label+'_'+key+'_difference'] = x-y
            rows[key] = float(np.linalg.norm(x-y)/np.linalg.norm(y))
        controls[label] = rows
    quadrature = []
    for i, case in enumerate((base, fine, small)):
        if case['max_induction_residual'] > 1e-10: raise ValueError('Original strict bound violated')
        low, high = case['quadrature']; differences = []
        for group in ('rate_errors', 'initial_errors'):
            for name in low[group]:
                for region in ('relative_L2', 'axis_element', 'wall_element'):
                    a = low[group][name][region]; b = high[group][name][region]
                    differences.append(dict(group=group, name=name, region=region, low=a, high=b,
                        absolute=abs(a-b), relative=abs(a-b)/max(abs(a), abs(b), np.finfo(float).tiny)))
        quadrature.append(dict(index=i, largest_relative=max(differences, key=lambda row: row['relative']),
                               largest_absolute=max(differences, key=lambda row: row['absolute'])))
    identity = provenance(['examples.split_inverse_comparison'])
    destination.mkdir(parents=True); freeze(destination, identity)
    path = destination/'differences.npz'; np.savez_compressed(path, **arrays)
    result = dict(terminal=True, source_sha256=identity['source_sha256'], input_hashes=inputs,
                  baseline_input_and_full_fluid_bitwise_equal=True, coefficient_relative_differences=controls,
                  coefficient_metric='Euclidean coefficient norms, not physical reconstruction norms',
                  quadrature_controls=quadrature, arrays_file=path.name, arrays_sha256=sha(path))
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print('split-inverse controls', controls, flush=True)
