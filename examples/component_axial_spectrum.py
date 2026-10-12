"""Exact discrete cell-Fourier energy partition of complete physical errors."""
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
    if specs[0] != specs[1]: raise ValueError('Same physical mesh and settings required')
    if any(summaries[0][k] != summaries[1][k] for k in ('reference_directory', 'reference_index')):
        raise ValueError('Same independent continuum mode required')
    mode = ArchivedContactMode(summaries[0]['reference_directory'], summaries[0]['reference_index'])
    if any(c['input_sha256'] != mode.input_sha256 or c['input_summary_sha256'] != mode.input_summary_sha256 for c in cases):
        raise ValueError('Unchanged independent reference required')
    if any(not np.array_equal(raw[0][k], raw[1][k]) for k in ('state', 'fluid_perturbation', 'magnetic_perturbation')):
        raise ValueError('Same physical base and tangent required')
    spec = specs[0]; ref = reference(dict(spec, kind='contact'))
    d = RegularCoordinateFluid(ref.radius, ref.length, spec['elements'], spec['points'],
        ref.ion_mass+electron_mass, degree=spec.get('degree', 5), gammas=(ref.gamma_e, ref.gamma_i))
    mesh = MHDMesh(ref.radius, .5*ref.length, spec['field_nr'], spec['field_nz'])
    rates = [{sign: SimpleNamespace(fluid=a[sign+'_fluid_rate'], magnetic=a[sign+'_magnetic_rate'])
        for sign in ('minus', 'plus')} for a in raw]
    coefficients = [{name: (a['plus_'+name]-a['minus_'+name])/(2*c['epsilon'])
        for name in ('electric', 'current_rate')} for a, c in zip(raw, cases)]
    identity = provenance(['examples.component_axial_spectrum'])
    destination.mkdir(parents=True); freeze(destination, identity)
    result = dict(terminal=False, source_sha256=identity['source_sha256'],
        inputs=[list(first), list(second)], input_hashes={p: h for item in loaded for p, h in item[3].items()},
        amplitudes=[c['epsilon'] for c in cases], quadrature=[])
    def save(): (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    save(); arrays = {}; wave = np.rint(np.fft.fftfreq(mesh.nz)*mesh.nz).astype(int)
    groups = {'all': np.ones(mesh.nz, bool), 'zero': wave == 0,
        'reference_band': abs(wave) <= mode.points//2,
        'above_reference_band': abs(wave) > mode.points//2,
        'upper_quarter': abs(wave) >= 3*mesh.nz/8}
    names = (*NAMES, *(group+'_'+axis for group in ('electric', 'current') for axis in ('radial', 'azimuthal', 'axial')))
    for order in (16, 20):
        print('Full physical component cell-Fourier audit, quadrature', order, flush=True)
        sample = TransportSampler(d, mesh, order); prints = [FieldFingerprints(len(sample.r)) for _ in cases]
        terms = ('first_error', 'second_error', 'difference', 'reference')
        spectral = {name: {key: [] for key in terms} for name in names}
        direct = {name: {key: [] for key in terms} for name in names}
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
                values = dict(first_error=fields[0][name]-expected[name], second_error=fields[1][name]-expected[name],
                    difference=fields[0][name]-fields[1][name], reference=expected[name])
                for key, value in values.items():
                    transformed = np.fft.fft(value.reshape(len(s.r), mesh.nz, order), axis=1)/mesh.nz
                    spectrum = mesh.nz*np.sum(abs(transformed)**2*s.weight[:, None, :order], axis=2)
                    spectral[name][key].append(np.sum(spectrum, axis=0))
                    direct[name][key].extend(np.sum(value*value*s.weight, axis=1).tolist())
        fingerprints = [p.result() for p in prints]
        for case, observed in zip(cases, fingerprints):
            original = next(q for q in case['quadrature'] if q['order'] == order)['field_fingerprints']
            if any(original[k] != v for k, v in observed.items()): raise ValueError('Original full physical fingerprints differ')
        metrics = {}
        for name in names:
            metrics[name] = {}
            for key in terms:
                spectrum = np.array([math.fsum(col) for col in np.asarray(spectral[name][key]).T])
                total = math.fsum(direct[name][key]); parseval = math.fsum(spectrum)
                relative_defect = abs(parseval-total)/max(total, np.finfo(float).tiny)
                if relative_defect > 256*np.finfo(float).eps: raise ValueError('Discrete Parseval identity failed')
                metrics[name][key] = dict(direct_squared_integral=total, spectral_squared_integral=parseval,
                    parseval_relative_defect=relative_defect,
                    fractions={group: math.fsum(spectrum[mask])/max(total, np.finfo(float).tiny) for group, mask in groups.items()},
                    largest_signed_mode=int(wave[np.argmax(spectrum)]))
                arrays[f'q{order}_{name}_{key}_spectrum'] = spectrum
        result['quadrature'].append(dict(order=order, metrics=metrics, original_fingerprints_reproduced=True,
            reconstructed_fingerprints=fingerprints)); save()
    arrays['signed_cell_mode'] = wave
    path = destination/'spectra.npz'; np.savez_compressed(path, **arrays)
    if provenance(['examples.component_axial_spectrum']) != identity: raise ValueError('Source changed')
    result.update(terminal=True, arrays_file=path.name, arrays_sha256=sha(path),
        scope='Cell-translation Fourier modes with intra-cell Gauss coordinate retained; aliases separated by nz fold together. No state filtering or plasma solve.')
    save(); return result
