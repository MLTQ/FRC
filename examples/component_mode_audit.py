"""Physical component accuracy and amplitude sensitivity on fixed radial bands."""
import json
import math
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.constants import electron_mass
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from .regular_plasma_reference import reference
from .archived_contact_mode import ArchivedContactMode
from .current_amplitude_audit import load
from .current_transport_sampling import TransportSampler
from .continuum_wall_probe import physical
from .radial_term_sampling import radial_slice
from .streamed_rate_audit import NAMES, FieldFingerprints
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, first, second):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    loaded = [load(*item) for item in (first, second)]
    summaries, cases, raw = [[row[k] for row in loaded] for k in range(3)]
    specs = [{k: v for k, v in c['spec'].items() if k not in ('epsilon', 'tolerance')} for c in cases]
    if specs[0] != specs[1]: raise ValueError('Same grid and physical settings required')
    if any(summaries[0][k] != summaries[1][k] for k in ('reference_directory', 'reference_index')):
        raise ValueError('Same independent reference required')
    mode = ArchivedContactMode(summaries[0]['reference_directory'], summaries[0]['reference_index'])
    if any(c['input_sha256'] != mode.input_sha256 or c['input_summary_sha256'] != mode.input_summary_sha256 for c in cases):
        raise ValueError('Unchanged independent mode required')
    if any(not np.array_equal(raw[0][k], raw[1][k]) for k in ('state', 'fluid_perturbation', 'magnetic_perturbation')):
        raise ValueError('Identical physical base and tangent required')
    spec = specs[0]; ref = reference(dict(spec, kind='contact'))
    d = RegularCoordinateFluid(ref.radius, ref.length, spec['elements'], spec['points'],
        ref.ion_mass+electron_mass, degree=spec.get('degree', 5), gammas=(ref.gamma_e, ref.gamma_i))
    mesh = MHDMesh(ref.radius, .5*ref.length, spec['field_nr'], spec['field_nz'])
    rates = [{sign: SimpleNamespace(fluid=a[sign+'_fluid_rate'], magnetic=a[sign+'_magnetic_rate'])
        for sign in ('minus', 'plus')} for a in raw]
    coefficients = [{name: (a['plus_'+name]-a['minus_'+name])/(2*c['epsilon'])
        for name in ('electric', 'current_rate')} for a, c in zip(raw, cases)]
    identity = provenance(['examples.component_mode_audit'])
    destination.mkdir(parents=True); freeze(destination, identity)
    names = (*NAMES, 'electric', 'current_rate', *(group+'_'+axis for group in ('electric', 'current') for axis in ('radial', 'azimuthal', 'axial')))
    arrays = {}; records = []
    for order in (16, 20):
        print('Complete component reconstruction quadrature', order, flush=True)
        sample = TransportSampler(d, mesh, order)
        prints = [FieldFingerprints(len(sample.r)) for _ in cases]
        rows = {name: {key: [] for key in ('first_error', 'second_error', 'difference', 'reference')} for name in names}
        maxima = {name: dict.fromkeys(rows[name], 0.) for name in names}
        for start in range(0, len(sample.r), 2*order):
            s = radial_slice(sample, start, start+2*order); target = mode.fields(s.R, s.Z)
            expected = dict(physical_rate=mode.growth*target['physical'], electric=target['electric'],
                current_rate=mode.growth*target['current'])
            fields = []
            for i, (a, case, rate, coefficient) in enumerate(zip(raw, cases, rates, coefficients)):
                values = {name: s.electric(v) for name, v in coefficient.items()}
                values['physical_rate'] = (physical(s, a['plus_state'], rate['plus'])-
                    physical(s, a['minus_state'], rate['minus']))/(2*case['epsilon'])
                prints[i].add(**values, target_physical=expected['physical_rate'],
                    target_electric=expected['electric'], target_current_rate=expected['current_rate'], weight=s.weight)
                fields.append(values)
            for values in (expected, *fields):
                values.update({name: values['physical_rate'][..., k] for k, name in enumerate(NAMES)})
                for source, prefix in [('electric', 'electric'), ('current_rate', 'current')]:
                    values.update({prefix+'_'+axis: values[source][..., k]
                        for k, axis in enumerate(('radial', 'azimuthal', 'axial'))})
            for name in names:
                terms = dict(first_error=fields[0][name]-expected[name], second_error=fields[1][name]-expected[name],
                    difference=fields[0][name]-fields[1][name], reference=expected[name])
                for key, value in terms.items():
                    square = np.sum(value*value, axis=-1) if value.ndim == 3 else value*value
                    rows[name][key].extend(np.sum(s.weight*square, axis=1).tolist())
                    maxima[name][key] = max(maxima[name][key], float(np.sqrt(np.max(square))))
        fingerprints = [item.result() for item in prints]
        for case, observed in zip(cases, fingerprints):
            original = next(q for q in case['quadrature'] if q['order'] == order)['field_fingerprints']
            if any(original[k] != v for k, v in observed.items()):
                raise ValueError('Reconstruction differs from original physical fields')
        masks = dict(volume=np.ones(len(sample.r), bool), axis=sample.r < .05*ref.radius,
            wall=sample.r >= .95*ref.radius, interior=(sample.r >= .05*ref.radius)&(sample.r < .95*ref.radius))
        metrics = {}
        for name, terms in rows.items():
            terms = {k: np.asarray(v) for k, v in terms.items()}
            totals = {key: math.fsum(value) for key, value in terms.items()}
            regions = {}
            for region, mask in masks.items():
                integral = {key: math.fsum(value[mask]) for key, value in terms.items()}
                regions[region] = dict(squared_integrals=integral,
                    relative_L2={key: float(np.sqrt(integral[key]/integral['reference'])) if integral['reference'] > 0 else None
                        for key in ('first_error', 'second_error', 'difference')},
                    fraction_of_domain_squared={key: integral[key]/totals[key] if totals[key] > 0 else None for key in totals})
            metrics[name] = dict(regions=regions, sampled_maxima=maxima[name])
            arrays.update({f'q{order}_{name}_{key}_radial_integrals': value for key, value in terms.items()})
        arrays[f'q{order}_radius'] = sample.r
        records.append(dict(order=order, metrics=metrics, reconstructed_fingerprints=fingerprints,
            original_fingerprints_reproduced=True))
    path = destination/'profiles.npz'; np.savez_compressed(path, **arrays)
    if provenance(['examples.component_mode_audit']) != identity: raise ValueError('Source changed')
    result = dict(terminal=True, source_sha256=identity['source_sha256'], inputs=[list(first), list(second)],
        amplitudes=[c['epsilon'] for c in cases], input_hashes={p: h for item in loaded for p, h in item[3].items()},
        regions=dict(axis='r < 0.05 R', wall='r >= 0.95 R', interior='0.05 R <= r < 0.95 R'),
        quadrature=records, arrays_file=path.name, arrays_sha256=sha(path),
        scope='Complete centered physical component rates and their difference; no plasma solve or acceptance change')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
