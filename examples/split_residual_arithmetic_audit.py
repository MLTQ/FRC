"""Independent exact electric-sum and vector attribution on captured iterates."""
from decimal import Decimal, localcontext
import json
from pathlib import Path
import numpy as np
from .implicit_plasma_records import provenance, freeze, sha


TERMS = ('metric_association', 'inertia_paths', 'inertia_curl_distribution',
    'inertia_dual_distribution', 'electric_curl_distribution',
    'electric_dual_distribution', 'rate_dual_distribution')


def run(destination, original='outputs/split_residual_probe'):
    root = Path(original); destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    record = json.loads((root/'summary.json').read_text())
    source = json.loads((root/'provenance.json').read_text())
    if not record['terminal'] or not record['observation_consistent'] or not record['matches_original_first_attempt_failure']:
        raise ValueError('Complete matched-rejection observation required')
    if record['source_sha256'] != source['source_sha256'] or any(sha(p) != h for p, h in source['source_files'].items()):
        raise ValueError('Unchanged recorded observation source required')
    path = root/record['arrays_file']
    if sha(path) != record['arrays_sha256']: raise ValueError('Unchanged full observation arrays required')
    data = dict(np.load(path)); identity = provenance(['examples.split_residual_arithmetic_audit'])
    destination.mkdir(parents=True); freeze(destination, identity)
    rows = []; arrays = {}; norm = np.linalg.norm
    for index in range(len(record['trace'])):
        prefix = f'refinement_{index}_'
        a, b, rounded, low = [data[prefix+key] for key in
            ('forcing_electric', 'inertia_electric', 'electric', 'electric_addition_low')]
        if not all(np.isfinite(v).all() for v in (a, b, rounded, low)):
            raise ValueError('Finite binary64 actual electric terms required')
        print('Exact Decimal actual-field audit refinement', index, flush=True)
        expected = np.empty_like(low); mismatch = 0
        with localcontext() as context:
            context.prec = 2200
            for k, (x, y, z, observed) in enumerate(zip(a, b, rounded, low)):
                exact = Decimal.from_float(float(x))+Decimal.from_float(float(y))-Decimal.from_float(float(z))
                expected[k] = float(exact)
                mismatch += exact != Decimal.from_float(float(observed))
        arrays[prefix+'decimal_addition_low'] = expected
        gap = data[prefix+'opposing_gap']; rhs = data[prefix+'rhs']; scale = max(norm(rhs), np.finfo(float).tiny)
        if not np.array_equal(gap, data[prefix+'operator_residual']+data[prefix+'physical_dual_residual']):
            raise ValueError('Saved opposing residual vector identity failed')
        terms = {key: data[prefix+key] for key in TERMS}
        if not np.array_equal(sum(terms.values()), data[prefix+'predicted_gap']):
            raise ValueError('Saved decomposition vector identity failed')
        denominator = max(float(gap@gap), np.finfo(float).tiny)
        gapnorm = max(norm(gap), np.finfo(float).tiny)
        attribution = {key: dict(norm_over_rhs=float(norm(value)/scale),
            signed_projection_on_gap=float(value@gap/denominator),
            gap_minus_term_fraction=float(norm(gap-value)/gapnorm)) for key, value in terms.items()}
        addition = data[prefix+'electric_addition_prediction']
        rows.append(dict(refinement=index, coefficients_checked=len(a), decimal_precision=2200,
            exact_decimal_mismatches=int(mismatch), bitwise_low_equal=bool(np.array_equal(expected, low)),
            gap_over_rhs=float(norm(gap)/scale), terms=attribution,
            predicted_remainder_fraction=float(norm(gap-sum(terms.values()))/gapnorm),
            gap_minus_electric_addition_prediction_fraction=float(norm(gap-addition)/gapnorm)))
    path = destination/'decimal_arrays.npz'; np.savez_compressed(path, **arrays)
    verified = all(row['exact_decimal_mismatches'] == 0 and row['bitwise_low_equal'] for row in rows)
    if provenance(['examples.split_residual_arithmetic_audit']) != identity: raise ValueError('Source changed')
    result = dict(terminal=True, exact_electric_sum_verified=verified, source_sha256=identity['source_sha256'],
        input_directory=str(root), input_summary_sha256=sha(root/'summary.json'), input_arrays_sha256=record['arrays_sha256'],
        observations=rows, arrays_file=path.name, arrays_sha256=sha(path),
        full_plasma_rate_accepted=False, scope='Observation of rejected original iterates, not solver repair')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    if not verified: raise ValueError('Exact actual-field Decimal audit disagrees with recorded two-sum terms')
    return result
