"""Direction-separated original-continuum rate controls with compact raw data."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from plasma_bench.harmonic_magnetic import HarmonicMagnetic
from .matched_mode_fixture import fixture
from .continuum_wall_probe import physical
from .current_transport_sampling import TransportSampler
from .regular_plasma_norms import norm
from .implicit_plasma_records import provenance, freeze, sha


def fingerprint(value):
    value = np.ascontiguousarray(value)
    header = json.dumps(dict(shape=value.shape, dtype=value.dtype.str), sort_keys=True).encode()
    return hashlib.sha256(header+value.tobytes()).hexdigest()


def measure(mode, d, system, u, du, db, epsilon, rates, order):
    sample = TransportSampler(d, system.geometry.mesh, order)
    target = mode.fields(sample.R, sample.Z)
    minus = u-epsilon*du; plus = u+epsilon*du
    differential = (physical(sample, plus, rates[1])-physical(sample, minus, rates[0]))/(2*epsilon)
    electric = sample.electric((rates[1].electric-rates[0].electric)/(2*epsilon))
    current = sample.electric((rates[1].current-rates[0].current)/(2*epsilon))
    names = ['density', 'radial_velocity', 'azimuthal_velocity', 'axial_velocity', 'electron_entropy', 'ion_entropy', 'radial_B', 'azimuthal_B', 'axial_B']
    expected = mode.growth*target['physical']
    errors = {name: norm(sample, differential[..., k], expected[..., k]) for k, name in enumerate(names)}
    errors['electric'] = norm(sample, electric, target['electric'])
    errors['current_rate'] = norm(sample, current, mode.growth*target['current'])
    n, _, _, dn, dv = sample.fluid(u, du)
    entropy = np.stack([sample.scalar(u[..., k])/n for k in (4, 5)], axis=-1)
    ds = np.stack([(sample.scalar(du[..., k])-entropy[..., k-4]*dn)/n for k in (4, 5)], axis=-1)
    initial = np.concatenate((dn[..., None], dv, ds, sample.magnetic(db)), axis=-1)
    initial_errors = {name: norm(sample, initial[..., k], target['physical'][..., k]) for k, name in enumerate(names)}
    fields = dict(physical_rate=differential, electric=electric, current_rate=current,
        target_physical=expected, target_electric=target['electric'], target_current_rate=mode.growth*target['current'],
        weight=sample.weight, initial_physical=initial)
    return dict(order=order, rate_errors=errors, initial_errors=initial_errors,
                field_fingerprints={k: fingerprint(v) for k, v in fields.items()})


def calculate(spec, epsilon=1e-5, orders=(16, 20)):
    start = time.monotonic(); mode, d, u, system, b, du, db = fixture(spec); rates = []
    if not np.isfinite(epsilon) or epsilon <= 0: raise ValueError('Positive perturbation amplitude required')
    arrays = dict(state=u, fluid_perturbation=du, magnetic_perturbation=db)
    for name, sign in [('minus', -1), ('plus', 1)]:
        state = u+sign*epsilon*du; psi, g = system.geometry.unpack(sign*epsilon*db)
        field = HarmonicMagnetic(system.geometry.mesh, psi, g, b.axial_bias)
        rate = system.derivative(state, field); rates.append(rate)
        arrays.update({name+'_'+key: value for key, value in dict(state=state, magnetic=field.packed_increment(system.geometry),
            fluid_rate=rate.fluid, magnetic_rate=rate.magnetic, electric=rate.electric, current_rate=rate.current).items()})
    measured = [measure(mode, d, system, u, du, db, epsilon, rates, order) for order in orders]
    result = dict(spec=spec, epsilon=epsilon, input_sha256=mode.input_sha256, continuum_growth=mode.growth,
        quadrature=measured, max_induction_residual=max(r.residual for r in rates),
        rate_diagnostics=[r.diagnostics for r in rates], elapsed_seconds=time.monotonic()-start)
    return result, arrays


def run(destination, specs):
    destination = Path(destination)
    if destination.exists(): raise ValueError('Refusing overwrite')
    seeds = ['examples.matched_mode_rates']; identity = provenance(seeds)
    destination.mkdir(parents=True); freeze(destination, identity); results = []
    for spec in specs:
        epsilon = spec.get('epsilon', 1e-5)
        result, arrays = calculate(spec, epsilon)
        raw = destination/f'case_{len(results):02}.npz'; np.savez_compressed(raw, **arrays)
        result.update(arrays_file=raw.name, arrays_sha256=sha(raw)); results.append(result)
        (destination/f'case_{len(results)-1:02}.json').write_text(json.dumps(result, indent=2)+'\n')
        errors = result['quadrature'][0]['rate_errors']
        print('matched rates', spec, 'Br', errors['radial_B']['relative_L2'], 'Jt', errors['current_rate']['relative_L2'], flush=True)
    if provenance(seeds) != identity: raise RuntimeError('Source changed')
    (destination/'summary.json').write_text(json.dumps(dict(terminal=True, source_sha256=identity['source_sha256'], cases=results), indent=2)+'\n')
