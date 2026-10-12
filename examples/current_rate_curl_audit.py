"""Separate magnetic-rate error from weak-versus-strong current reconstruction."""
import json
import math
from pathlib import Path
import numpy as np
from scipy.constants import electron_mass, mu_0
from plasma_bench.regular_coordinate_fluid import RegularCoordinateFluid
from plasma_bench.mhd_mesh import MHDMesh
from .regular_plasma_reference import reference
from .archived_contact_mode import ArchivedContactMode
from .current_amplitude_audit import load
from .current_transport_sampling import TransportSampler, bases
from .radial_term_sampling import radial_slice
from .streamed_rate_audit import FieldFingerprints
from .implicit_plasma_records import provenance, freeze, sha


def run(destination, cases):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    loaded = [load(*item) for item in cases]
    identity = provenance(['examples.current_rate_curl_audit'])
    destination.mkdir(parents=True); freeze(destination, identity)
    result = dict(terminal=False, source_sha256=identity['source_sha256'],
        inputs=[list(item) for item in cases], input_hashes={p: h for row in loaded for p, h in row[3].items()}, cases=[])
    def save(): (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    save(); arrays = {}
    for index, (summary, case, raw, _) in enumerate(loaded):
        spec = case['spec']; ref = reference(dict(spec, kind='contact'))
        mode = ArchivedContactMode(summary['reference_directory'], summary['reference_index'])
        if mode.input_sha256 != case['input_sha256'] or mode.input_summary_sha256 != case['input_summary_sha256']:
            raise ValueError('Unchanged independent mode required')
        d = RegularCoordinateFluid(ref.radius, ref.length, spec['elements'], spec['points'], ref.ion_mass+electron_mass,
            degree=spec.get('degree', 5), gammas=(ref.gamma_e, ref.gamma_i))
        mesh = MHDMesh(ref.radius, .5*ref.length, spec['field_nr'], spec['field_nz'])
        magnetic = (raw['plus_magnetic_rate']-raw['minus_magnetic_rate'])/(2*case['epsilon'])
        current = (raw['plus_current_rate']-raw['minus_current_rate'])/(2*case['epsilon'])
        nodes = (mesh.nr+1)*mesh.nz
        psi = magnetic[:nodes].reshape(mesh.nr+1, mesh.nz); g = magnetic[nodes:].reshape(mesh.nr, mesh.nz)
        records = []
        for order in (16, 20):
            print('Independent full magnetic/current-rate curl comparison', index, order, flush=True)
            sample = TransportSampler(d, mesh, order); fingerprint = FieldFingerprints(len(sample.r))
            dzz = bases(mesh.zf, sample.z, 5, True)[2]
            terms = ('current_error', 'strong_error', 'weak_minus_strong', 'reference', 'cross')
            integrals = {axis: {key: [] for key in terms} for axis in ('radial', 'azimuthal', 'axial')}
            for start in range(0, len(sample.r), 2*order):
                s = radial_slice(sample, start, start+2*order)
                j = s.electric(current); target = mode.growth*mode.fields(s.R, s.Z)['current']
                dxx = bases(.5*mesh.rf**2, .5*s.r**2, 5)[2]
                strong = np.stack((-s.R*(s.ex@g@s.dez.T),
                    -(s.nx@psi@dzz.T)/s.R-s.R*(dxx@psi@s.nz.T),
                    (2*s.ex+s.R**2*s.dex)@g@s.ez.T), axis=-1)/mu_0
                fingerprint.add(current_rate=j, target_current_rate=target, weight=s.weight)
                first = j-target; second = strong-target; defect = j-strong
                for k, axis in enumerate(integrals):
                    values = dict(current_error=first[..., k]**2, strong_error=second[..., k]**2,
                        weak_minus_strong=defect[..., k]**2, reference=target[..., k]**2,
                        cross=2*second[..., k]*defect[..., k])
                    for key, value in values.items():
                        integrals[axis][key].extend(np.sum(s.weight*value, axis=1).tolist())
            observed = fingerprint.result()
            original = next(q for q in case['quadrature'] if q['order'] == order)['field_fingerprints']
            if any(original[k] != v for k, v in observed.items()): raise ValueError('Original complete current fields differ')
            metrics = {}
            for axis, values in integrals.items():
                totals = {key: math.fsum(value) for key, value in values.items()}
                decomposition = totals['strong_error']+totals['weak_minus_strong']+totals['cross']
                metrics[axis] = dict(squared_integrals=totals,
                    relative_L2={key: float(np.sqrt(totals[key]/totals['reference']))
                        for key in ('current_error', 'strong_error', 'weak_minus_strong')},
                    energy_decomposition_relative_defect=abs(totals['current_error']-decomposition)/max(totals['current_error'],np.finfo(float).tiny))
                arrays.update({f'case{index}_q{order}_{axis}_{key}_rows': np.asarray(value) for key, value in values.items()})
            records.append(dict(order=order, metrics=metrics, original_fingerprints_reproduced=True,
                reconstructed_fingerprints=observed))
        result['cases'].append(dict(spec=spec, epsilon=case['epsilon'], quadrature=records)); save()
    path = destination/'curl_comparison.npz'; np.savez_compressed(path, **arrays)
    if provenance(['examples.current_rate_curl_audit']) != identity: raise ValueError('Source changed')
    result.update(terminal=True, arrays_file=path.name, arrays_sha256=sha(path),
        scope='Strong curl inside polynomial cells excludes distributional face jumps; weak-minus-strong is not itself a failed discrete Ampere equation.')
    save(); return result
