"""Physical amplitude sensitivity of complete saved current/electric rates."""
import json
import math
from pathlib import Path
import numpy as np
from scipy.constants import electron_mass
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from .regular_plasma_reference import reference
from .archived_contact_mode import ArchivedContactMode
from .current_transport_sampling import TransportSampler
from .radial_term_sampling import radial_slice
from .streamed_rate_audit import FieldFingerprints
from .implicit_plasma_records import provenance, freeze, sha


def load(directory, index):
    directory = Path(directory)
    summary = json.loads((directory/'summary.json').read_text())
    source = json.loads((directory/'provenance.json').read_text())
    if not summary['terminal'] or summary['source_sha256'] != source['source_sha256']:
        raise ValueError('Terminal original source-identified workflow required')
    if any(sha(p) != h for p, h in source['source_files'].items()):
        raise ValueError('Unchanged original calculation source required')
    case = summary['cases'][index]
    if case['status'] != 'complete':
        raise ValueError('A failed or unfinished solve has no physical rate')
    if not np.isfinite(case['max_induction_residual']) or case['max_induction_residual'] > case['spec'].get('tolerance', 1e-10):
        raise ValueError('Original declared strict bound required')
    path = directory/case['arrays_file']
    if sha(path) != case['arrays_sha256']:
        raise ValueError('Unchanged raw arrays required')
    return summary, case, dict(np.load(path)), {str(path): sha(path),
        str(directory/'summary.json'): sha(directory/'summary.json'),
        str(directory/'provenance.json'): sha(directory/'provenance.json')}


def run(destination, first, second):
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Refusing overwrite')
    loaded = [load(*item) for item in (first, second)]
    summaries, cases, raw = [[row[k] for row in loaded] for k in range(3)]
    spec = {k: v for k, v in cases[0]['spec'].items() if k not in ('epsilon', 'tolerance')}
    if spec != {k: v for k, v in cases[1]['spec'].items() if k not in ('epsilon', 'tolerance')}:
        raise ValueError('Same physical grid and settings required')
    for key in ('reference_directory', 'reference_index'):
        if summaries[0][key] != summaries[1][key]:
            raise ValueError('Same independent reference required')
    mode = ArchivedContactMode(summaries[0]['reference_directory'], summaries[0]['reference_index'])
    for case in cases:
        if case['input_sha256'] != mode.input_sha256 or case['input_summary_sha256'] != mode.input_summary_sha256:
            raise ValueError('Unchanged independent reference required')
    for name in ('state', 'fluid_perturbation', 'magnetic_perturbation'):
        if not np.array_equal(raw[0][name], raw[1][name]):
            raise ValueError('Same base state and physical perturbation required')
    ref = reference(dict(spec, kind='contact'))
    d = RegularCoordinateFluid(ref.radius, ref.length, spec['elements'], spec['points'],
        ref.ion_mass+electron_mass, degree=spec.get('degree', 5), gammas=(ref.gamma_e, ref.gamma_i))
    mesh = MHDMesh(ref.radius, .5*ref.length, spec['field_nr'], spec['field_nz'])
    coefficients = [{name: (a['plus_'+key]-a['minus_'+key])/(2*case['epsilon'])
        for name, key in [('current_rate', 'current_rate'), ('electric', 'electric')]}
        for a, case in zip(raw, cases)]
    identity = provenance(['examples.current_amplitude_audit'])
    destination.mkdir(parents=True); freeze(destination, identity)
    records = []; arrays = {}
    for order in (16, 20):
        print('physical amplitude difference quadrature', order, flush=True)
        sample = TransportSampler(d, mesh, order)
        prints = [FieldFingerprints(len(sample.r)) for _ in cases]
        names = ('current_rate', 'current_radial', 'current_azimuthal', 'current_axial', 'electric')
        rows = {name: {key: [] for key in ('difference', 'reference', 'first', 'second')} for name in names}
        maxima = dict.fromkeys(names, 0.)
        for start in range(0, len(sample.r), 2*order):
            s = radial_slice(sample, start, start+2*order)
            target = mode.fields(s.R, s.Z)
            expected = dict(current_rate=mode.growth*target['current'], electric=target['electric'])
            fields = [{name: s.electric(value) for name, value in c.items()} for c in coefficients]
            for fingerprint, values in zip(prints, fields):
                fingerprint.add(**values, target_current_rate=expected['current_rate'],
                    target_electric=expected['electric'], weight=s.weight)
            for component, name in enumerate(names[1:4]):
                expected[name] = expected['current_rate'][..., component]
                for values in fields:
                    values[name] = values['current_rate'][..., component]
            for name in names:
                values = dict(difference=fields[0][name]-fields[1][name],
                    reference=expected[name], first=fields[0][name], second=fields[1][name])
                for key, value in values.items():
                    square = np.sum(value*value, axis=-1) if value.ndim == 3 else value*value
                    rows[name][key].extend(np.sum(s.weight*square, axis=1).tolist())
                    if key == 'difference':
                        maxima[name] = max(maxima[name], float(np.sqrt(np.max(square))))
        fingerprints = [item.result() for item in prints]
        for case, observed in zip(cases, fingerprints):
            original = next(q for q in case['quadrature'] if q['order'] == order)['field_fingerprints']
            if any(original[name] != digest for name, digest in observed.items()):
                raise ValueError('Reconstruction differs from original physical audit')
        metrics = {}
        for name, terms in rows.items():
            integral = {key: math.fsum(value) for key, value in terms.items()}
            metrics[name] = dict(max_absolute_difference=maxima[name], squared_integrals=integral,
                **{'relative_to_'+key: float(np.sqrt(integral['difference']/max(integral[key], np.finfo(float).tiny)))
                   for key in ('reference', 'first', 'second')})
            arrays.update({f'q{order}_{name}_{key}_radial_integrals': np.asarray(value) for key, value in terms.items()})
        records.append(dict(order=order, metrics=metrics, original_fingerprints_reproduced=True,
            reconstructed_fingerprints=fingerprints))
    path = destination/'differences.npz'; np.savez_compressed(path, **arrays)
    if provenance(['examples.current_amplitude_audit']) != identity:
        raise ValueError('Source changed during audit')
    result = dict(terminal=True, source_sha256=identity['source_sha256'], inputs=[list(first), list(second)],
        input_hashes={p: h for item in loaded for p, h in item[3].items()},
        amplitudes=[case['epsilon'] for case in cases], quadrature=records,
        quadrature_absolute_differences={name: abs(records[0]['metrics'][name]['relative_to_reference']-
            records[1]['metrics'][name]['relative_to_reference']) for name in names},
        arrays_file=path.name, arrays_sha256=sha(path),
        scope='Physical volume current-rate and electric differences, not coefficient norms or error-norm differences')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
