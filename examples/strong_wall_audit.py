"""Independent physical quadrature and optional complete immutable replay."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.regular_plasma_step import field_from_packed
from .current_transport_sampling import TransportSampler
from .strong_wall_reference import initialize
from .regular_plasma_reference import exact
from .regular_plasma_evolution import field_norms
from .regular_plasma_norms import norm
from .implicit_plasma_records import provenance, freeze, sha
from .strong_wall_evolution import run_case


def measure(spec, arrays, time, order):
    ref, d, _, system, initial_field = initialize(spec)
    sample = TransportSampler(d, system.geometry.mesh, order)
    values = []; norms = {}
    for prefix, when in (('initial', 0.), ('final', time)):
        u = arrays[prefix+'_state']; packed = arrays[prefix+'_magnetic']
        field = field_from_packed(system.geometry, packed, initial_field.axial_bias)
        n, v, _ = sample.fluid(u)
        entropy = np.stack([sample.scalar(u[..., k])/n for k in (4, 5)], axis=-1)
        h = E*np.exp(entropy)*(n[..., None]/d.thermal.n0)**(d.thermal.gammas-1)/(d.thermal.gammas-1)
        model = RegularCurrent(system.mapping, u[..., 0], ref.ion_mass)
        current = sample.electric(field.current(model)); magnetic = sample.magnetic(packed)
        values.append((n, v, entropy, h, current, magnetic))
        norms[prefix] = field_norms(spec, ref, d, u, system, field, when, order)
        rate = system.derivative(u, field)
        target = exact(ref, spec['kind'], sample.R, sample.Z, when)
        norms[prefix]['electric'] = norm(sample, sample.electric(rate.electric), target['electric'])
        norms[prefix]['induction_residual'] = rate.residual
        norms[prefix]['wall_trace_relative'] = system.geometry.wall_residual(field.g)
        norms[prefix]['unprojected_curl_relative'] = rate.diagnostics['induction_residuals']['unprojected_curl_relative']
    n0, v0, s0, h0, j0, b0 = values[0]
    n, v, s, h, j, b = values[1]
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
    full_b0 = b0.copy(); full_b0[..., 2] += initial_field.axial_bias
    initial_energy = float(np.sum(sample.weight*(.5*d.mass*n0*np.sum(v0*v0, axis=-1)
        +n0*h0.sum(axis=-1)+.5*model.reduced_mass/E**2*np.sum(j0*j0, axis=-1)/n0
        +np.sum(full_b0*full_b0, axis=-1)/(2*MU0))))
    total = sum(changes.values()); particles = float(np.sum(sample.weight*n0))
    result = dict(order=order, norms=norms, energy_changes=changes, total_energy_change=total,
        initial_energy=initial_energy, relative_energy_drift=total/initial_energy,
        particle_relative_drift=float(np.sum(sample.weight*dn)/particles),
        entropy_change_per_particle=[float(np.sum(sample.weight*dc[..., k])/particles) for k in range(2)],
        point_density_min=float(arrays['final_state'][..., 0].min()), sampled_density_min=float(n.min()),
        temperature_min=np.min(h*(d.thermal.gammas-1)/E, axis=(0, 1)).tolist(),
        temperature_max=np.max(h*(d.thermal.gammas-1)/E, axis=(0, 1)).tolist())
    if spec['kind'] == 'wave': result['change_over_wave_energy'] = total/ref.perturbation_energy()
    return result


def without_runtime(value):
    if isinstance(value, dict): return {k: without_runtime(v) for k, v in value.items() if k != 'elapsed_seconds'}
    if isinstance(value, list): return [without_runtime(v) for v in value]
    return value


def run(directory, destination, replay=False):
    directory = Path(directory); destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing to overwrite independent audit')
    original = json.loads((directory/'provenance.json').read_text())
    root = Path(__file__).resolve().parents[1]
    if any(sha(root/name) != digest for name, digest in original['source_files'].items()):
        raise ValueError('Original calculation source required for audit/replay')
    seeds = ['examples.strong_wall_audit', 'examples.strong_wall_run']
    identity = provenance(seeds)
    if any(identity[k] != original[k] for k in ('python', 'numpy', 'scipy')):
        raise ValueError('Original runtime required')
    record = json.loads((directory/'case.json').read_text())
    raw = directory/record['arrays_file']
    if sha(raw) != record['arrays_sha256']: raise ValueError('Original arrays hash required')
    arrays = dict(np.load(raw)); destination.mkdir(parents=True); freeze(destination, identity)
    result = dict(input_directory=str(directory), input_source_sha256=original['source_sha256'],
        arrays_sha256=sha(raw), source_sha256=identity['source_sha256'], spec=record['spec'],
        order16=measure(record['spec'], arrays, record['duration'], 16),
        order20=measure(record['spec'], arrays, record['duration'], 20), full_replay=False)
    if replay:
        repeated, repeated_arrays = run_case(record['spec'])
        same = {k: np.array_equal(v, repeated_arrays[k]) for k, v in arrays.items()}
        expected = {k: record[k] for k in repeated}
        metrics = without_runtime(repeated) == without_runtime(expected)
        result.update(full_replay=True, bitwise_arrays=same, bitwise_metrics=metrics)
        if not all(same.values()) or not metrics: raise RuntimeError('Complete trajectory replay differed')
    if provenance(seeds) != identity: raise RuntimeError('Audit source changed')
    result['terminal'] = True
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print('independent audit complete', directory, 'replay', replay, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('directory'); parser.add_argument('destination')
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args(); run(args.directory, args.destination, args.replay)
