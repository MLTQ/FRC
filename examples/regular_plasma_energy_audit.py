"""Separate continuous physical energy audit of frozen coupled trajectories."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.regular_plasma_step import field_from_packed
from .current_transport_sampling import TransportSampler
from .regular_plasma_reference import initialize
from .regular_plasma_records import provenance, freeze, sha


def measure(spec, arrays, order):
    ref, d, _, system, initial_field = initialize(spec)
    sample = TransportSampler(d, system.geometry.mesh, order)
    values = []
    for prefix in ('initial', 'final'):
        u = arrays[prefix+'_state']
        field = field_from_packed(system.geometry, arrays[prefix+'_magnetic'], initial_field.axial_bias)
        n, v, _ = sample.fluid(u)
        entropy = np.stack([sample.scalar(u[..., k])/n for k in (4, 5)], axis=-1)
        h = E*np.exp(entropy)*(n[..., None]/d.thermal.n0)**(d.thermal.gammas-1)/(d.thermal.gammas-1)
        model = RegularCurrent(system.mapping, u[..., 0], ref.ion_mass)
        current = sample.electric(field.current(model))
        magnetic = sample.magnetic(arrays[prefix+'_magnetic'])
        values.append((n, v, entropy, h, current, magnetic))
    n0, v0, s0, h0, j0, b0 = values[0]
    n, v, s, h, j, b = values[1]
    # Apply the linear reconstruction to changes BEFORE subtracting large
    # reconstructed backgrounds. These are the same physical EOS differences.
    change = arrays['final_state']-arrays['initial_state']
    dn = sample.scalar(change[..., 0])
    dc = np.stack([sample.scalar(change[..., k]) for k in (4, 5)], axis=-1)
    ds = (dc-s0*dn[..., None])/n[..., None]
    thermal = n0[..., None]*h0*np.expm1(ds+d.thermal.gammas*np.log1p(dn/n0)[..., None])
    kinetic = .5*d.mass*(n*np.sum(v*v, axis=-1)-n0*np.sum(v0*v0, axis=-1))
    current = .5*model.reduced_mass/E**2*(np.sum(j*j, axis=-1)/n-np.sum(j0*j0, axis=-1)/n0)
    magnetic = (np.sum(b*b, axis=-1)-np.sum(b0*b0, axis=-1))/(2*MU0)
    cross = initial_field.axial_bias*sample.magnetic(arrays['final_magnetic']-arrays['initial_magnetic'])[..., 2]/MU0
    changes = {name: float(np.sum(sample.weight*value)) for name, value in
               [('bulk_kinetic', kinetic), ('electron_internal', thermal[..., 0]),
                ('ion_internal', thermal[..., 1]), ('current', current),
                ('magnetic_increment', magnetic), ('magnetic_bias_cross', cross)]}
    total = sum(changes.values())
    return dict(order=order, energy_changes=changes, total_change=total,
        change_over_wave_energy=total/ref.perturbation_energy(),
        particle_relative_drift=float(np.sum(sample.weight*dn)/np.sum(sample.weight*n0)),
        species_entropy_change=[float(np.sum(sample.weight*dc[..., k])) for k in range(2)],
        temperature_min=np.min(h*(d.thermal.gammas-1)/E, axis=(0, 1)).tolist(),
        temperature_max=np.max(h*(d.thermal.gammas-1)/E, axis=(0, 1)).tolist())


def run(directory, destination):
    directory = Path(directory); destination = Path(destination)
    if destination.exists():
        raise ValueError('Refusing to overwrite independent energy audit')
    identity = provenance()
    original = json.loads((directory/'provenance.json').read_text())
    if original != identity:
        raise ValueError('Original coupled source/runtime required')
    identity['source_files']['examples/regular_plasma_energy_audit.py'] = sha(__file__)
    identity['source_sha256'] = hashlib.sha256(json.dumps(identity['source_files'], sort_keys=True).encode()).hexdigest()
    destination.mkdir(parents=True); freeze(destination, identity)
    result = []
    for path in sorted(directory.glob('case_*.json')):
        record = json.loads(path.read_text())
        if record['kind'] != 'evolution':
            continue
        raw = directory/record['arrays_file']; assert sha(raw) == record['arrays_sha256']
        arrays = dict(np.load(raw))
        result.append(dict(file=path.name, arrays_sha256=sha(raw), spec=record['spec'],
                           grid_change=record['grid_change'],
                           order16=measure(record['spec'], arrays, 16), order20=measure(record['spec'], arrays, 20)))
    root = Path(__file__).resolve().parents[1]
    assert all(sha(root/name) == value for name, value in identity['source_files'].items())
    (destination/'summary.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],
        coupled_source_sha256=original['source_sha256'], terminal=True, cases=result), indent=2)+'\n')
    print(len(result), 'independent physical energy audits complete', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('directory'); parser.add_argument('destination')
    args = parser.parse_args(); run(args.directory, args.destination)
