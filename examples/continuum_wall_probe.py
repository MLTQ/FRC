"""Compare actual wall-RHS directional rates with an independent growing mode."""
import json
from pathlib import Path
import numpy as np
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from .strong_wall_reference import initialize
from .continuum_contact_mode import ContinuumContactMode
from .finite_mass_manufactured import cell_average
from .current_transport_sampling import TransportSampler
from .regular_plasma_norms import norm
from .implicit_plasma_records import provenance, freeze, sha


def fixture(elements):
    mode = ContinuumContactMode()
    spec = dict(kind='contact', elements=elements)
    _, d, u, system, b = initialize(spec); u[..., 1:4] = 0.
    target = mode.fields(d.r, d.z)['physical']; du = np.zeros_like(u); du[..., 0] = target[..., 0]
    np.divide(u[..., :1]*target[..., 1:3], d.r[..., None], out=du[..., 1:3], where=d.r[..., None] > 0)
    du[-1, -1, :, 1] = 0.; du[..., 3] = u[..., 0]*target[..., 3]
    du[..., 4:] = u[..., :1]*target[..., 4:6]+u[..., 4:]/u[..., :1]*du[..., :1]
    mesh = system.geometry.mesh
    psi = mode.fields(*mesh.corners)['psi']; psi[[0, -1]] = 0.; psi[:, -1] = psi[:, 0]
    g = cell_average(mesh, lambda r, z: mode.fields(r, z)['physical'][..., 7]/r, order=16)
    g = system.geometry.radial_projection@g
    db = system.geometry.pack(psi, g)
    return mode, d, u, system, b, du, db


def physical(sample, state, rate):
    n, _, _, nt, vt = sample.fluid(state, rate.fluid)
    entropy = np.stack([sample.scalar(state[..., k])/n for k in (4, 5)], axis=-1)
    st = np.stack([(sample.scalar(rate.fluid[..., k])-entropy[..., k-4]*nt)/n for k in (4, 5)], axis=-1)
    return np.concatenate((nt[..., None], vt, st, sample.magnetic(rate.magnetic)), axis=-1)


def calculate(elements, epsilon, order):
    mode, d, u, system, b, du, db = fixture(elements)
    sample = TransportSampler(d, system.geometry.mesh, order); target = mode.fields(sample.R, sample.Z)
    rates = []; physical_rates = []; electric = []; current = []
    for sign in (-1, 1):
        state = u+sign*epsilon*du
        psi, g = system.geometry.unpack(sign*epsilon*db)
        field = HarmonicMagnetic(system.geometry.mesh, psi, g, b.axial_bias)
        rate = system.derivative(state, field); rates.append(rate)
        physical_rates.append(physical(sample, state, rate)); electric.append(sample.electric(rate.electric))
        current.append(sample.electric(rate.current))
    derivative = (physical_rates[1]-physical_rates[0])/(2*epsilon)
    et = (electric[1]-electric[0])/(2*epsilon); jt = (current[1]-current[0])/(2*epsilon)
    names = ['density', 'radial_velocity', 'azimuthal_velocity', 'axial_velocity', 'electron_entropy', 'ion_entropy', 'radial_B', 'azimuthal_B', 'axial_B']
    expected = mode.growth*target['physical']
    errors = {name: norm(sample, derivative[..., k], expected[..., k]) for k, name in enumerate(names)}
    errors['electric'] = norm(sample, et, target['electric']); errors['current_rate'] = norm(sample, jt, mode.growth*target['current'])
    # Differential initial representation uses the independent quotient sampler.
    n, _, _, dn, dv = sample.fluid(u, du)
    entropy = np.stack([sample.scalar(u[..., k])/n for k in (4, 5)], axis=-1)
    ds = np.stack([(sample.scalar(du[..., k])-entropy[..., k-4]*dn)/n for k in (4, 5)], axis=-1)
    initial = np.concatenate((dn[..., None], dv, ds, sample.magnetic(db)), axis=-1)
    initial_errors = {name: norm(sample, initial[..., k], target['physical'][..., k]) for k, name in enumerate(names)}
    result = dict(elements=elements, epsilon=epsilon, order=order, input_sha256=mode.input_sha256,
        continuum_growth=mode.growth, rate_errors=errors, initial_errors=initial_errors,
        max_induction_residual=max(rate.residual for rate in rates))
    return result, dict(state=u, fluid_perturbation=du, magnetic_perturbation=db,
        physical_rate=derivative, electric=et, current_rate=jt, target_physical=expected,
        target_electric=target['electric'], target_current_rate=mode.growth*target['current'], weight=sample.weight)


def run(destination):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    seeds = ['examples.continuum_wall_probe']; identity = provenance(seeds)
    destination.mkdir(parents=True); freeze(destination, identity); results = []
    for elements, epsilon in [(1, 1e-5), (2, 1e-5), (4, 1e-5), (8, 1e-5), (16, 1e-5), (8, 1e-4)]:
        for order in (16, 20):
            result, arrays = calculate(elements, epsilon, order)
            raw = destination/f'case_{len(results):02}.npz'; np.savez_compressed(raw, **arrays)
            result.update(arrays_file=raw.name, arrays_sha256=sha(raw)); results.append(result)
            (destination/f'case_{len(results)-1:02}.json').write_text(json.dumps(result, indent=2)+'\n')
            print('continuum mode rates', elements, epsilon, order, result['rate_errors']['electric']['relative_L2'], flush=True)
    if provenance(seeds) != identity: raise RuntimeError('Source changed')
    (destination/'summary.json').write_text(json.dumps(dict(terminal=True, source_sha256=identity['source_sha256'], cases=results), indent=2)+'\n')


if __name__ == '__main__':
    import sys
    run(sys.argv[1])
