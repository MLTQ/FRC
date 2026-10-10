"""Fixed physical wall bands and error correlations from saved integrals."""
import json
from pathlib import Path
import numpy as np
from .implicit_plasma_records import provenance, freeze, sha


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    records = []; inputs = {}
    for directory in ('radial16_term_audit', 'radial16_drive_audit'):
        root = Path('outputs')/directory; path = root/'summary.json'
        summary = json.loads(path.read_text()); inputs[str(path)] = sha(path)
        if not summary['terminal']: raise ValueError('Terminal original audit required')
        for index, case in enumerate(summary['cases']):
            path = root/case['arrays_file']; inputs[str(path)] = sha(path)
            if sha(path) != case['arrays_sha256']: raise ValueError('Original integrals changed')
            raw = dict(np.load(path)); orders = []
            for order in (16, 20):
                prefix = f'q{order}_'; radius = raw[prefix+'radius']; count = len(radius)
                if count % 32: raise ValueError('Fixed physical band must align with cell faces')
                bands = {'axis_one_32': slice(0, count//32), 'interior': slice(count//32, 31*count//32),
                         'wall_one_32': slice(31*count//32, None)}
                terms = {}
                for key, value in raw.items():
                    if not key.startswith(prefix) or value.shape != (count, 3, 3): continue
                    error = value[:, 2]; total = error.sum(axis=0)
                    fractions = {}
                    for name, band in bands.items():
                        subset = error[band].sum(axis=0)
                        fractions[name] = [float(a/b) if b > 0 else None for a, b in zip(subset, total)]
                    for component in range(3):
                        if total[component] > 0 and not np.isclose(sum(v[component] for v in fractions.values()), 1, rtol=1e-13, atol=1e-13):
                            raise ValueError('Band partition does not reconstruct complete squared error')
                    terms[key[len(prefix):]] = fractions
                correlations = {}
                for group in ('momentum', 'faraday'):
                    key = prefix+group+'_error_cross'
                    if key not in raw: continue
                    matrix = raw[key].sum(axis=0); diagonal = np.diagonal(matrix, axis1=0, axis2=1).T
                    denominator = np.sqrt(diagonal[:, None, :]*diagonal[None, :, :])
                    correlations[group] = np.divide(matrix, denominator, out=np.zeros_like(matrix), where=denominator>0).tolist()
                orders.append(dict(order=order, squared_error_fractions=terms, error_correlations=correlations))
            records.append(dict(directory=directory, index=index, spec=case['spec'], quadrature=orders,
                                original_full_plasma_status=case.get('original_full_plasma_status', 'complete')))
    identity = provenance(['examples.radial_error_localization'])
    destination.mkdir(parents=True); freeze(destination, identity)
    result = dict(terminal=True, source_sha256=identity['source_sha256'], input_hashes=inputs,
                  fixed_wall_band='31/32 <= r/R <= 1 on every grid',
                  correlation_zero_denominator='reported as 0, not a measured nonzero correlation', cases=records)
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
