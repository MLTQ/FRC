"""Unforced full coupled trajectories against independent continuum fixtures."""
import time
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from plasma_bench.regular_current import RegularCurrent
from plasma_bench.regular_plasma_step import rk4
from .regular_plasma_reference import initialize, exact
from .current_transport_sampling import TransportSampler
from .regular_plasma_norms import norm


def field_norms(spec, ref, d, state, system, field, time_value, order):
    sample = TransportSampler(d, system.geometry.mesh, order)
    target = exact(ref, spec['kind'], sample.R, sample.Z, time_value)
    n, v, _ = sample.fluid(state)
    s = np.stack([sample.scalar(state[..., k])/n for k in (4, 5)], axis=-1)
    p = E*n[..., None]*np.exp(s)*(n[..., None]/d.thermal.n0)**(d.thermal.gammas-1)
    b = sample.magnetic(field.packed_increment(system.geometry))
    model = RegularCurrent(system.mapping, state[..., 0], ref.ion_mass)
    j = sample.electric(field.current(model))
    p0 = np.array([ref.electron_pressure, ref.ion_pressure]) if spec['kind'] == 'wave' else np.array([ref.electron_ev, ref.ion_ev])*ref.number_density*E
    norms = dict(density=norm(sample, n-ref.number_density, target['n']-ref.number_density),
        electron_pressure=norm(sample, p[..., 0]-p0[0], target['pressure'][..., 0]-p0[0]),
        ion_pressure=norm(sample, p[..., 1]-p0[1], target['pressure'][..., 1]-p0[1]))
    if spec['kind'] == 'wave':
        norms.update(velocity=norm(sample, v, target['velocity']),
                     magnetic=norm(sample, b, target['magnetic']), current=norm(sample, j, target['current']))
    else:
        norms.update(velocity=norm(sample, v, target['velocity']),
                     magnetic_max=float(np.max(np.linalg.norm(b, axis=-1))),
                     current_max=float(np.max(np.linalg.norm(j, axis=-1))))
    return norms


def grid_change(d, system, ref, initial, initial_field, state, field):
    f0, f = d.fields(initial), d.fields(state)
    n0, n = initial[..., 0], state[..., 0]
    thermal = n0[..., None]*f0['specific_internal']*np.expm1(
        f['entropy']-f0['entropy']+d.thermal.gammas*np.log1p((n-n0)/n0)[..., None])
    bulk = float(np.sum(d.volume*(f['kinetic']-f0['kinetic']+thermal.sum(axis=-1))))
    model0 = RegularCurrent(system.mapping, n0, ref.ion_mass)
    model = RegularCurrent(system.mapping, n, ref.ion_mass)
    current = float(np.sum(model.mesh.volumes*(model.current_energy(field.current(model))-
                                             model0.current_energy(initial_field.current(model0)))))
    magnetic = field.energy_components(model)['increment']-initial_field.energy_components(model0)['increment']
    return dict(fluid=bulk, current=current, magnetic=magnetic, total=bulk+current+magnetic,
                particles=float(np.sum(d.volume*(n-n0))),
                entropy=[float(np.sum(d.volume*(state[..., k]-initial[..., k]))) for k in (4, 5)],
                axial_momentum=float(d.mass*np.sum(d.volume*(state[..., 3]-initial[..., 3]))),
                angular_momentum=float(d.mass*np.sum(d.volume*2*d.x*(state[..., 2]-initial[..., 2]))))


def run_case(spec):
    start = time.monotonic()
    ref, d, initial, system, initial_field = initialize(spec)
    duration = ref.period if spec['kind'] == 'contact' else 2*np.pi/abs(ref.omega)
    duration *= spec.get('periods', 1.)
    count = spec.get('steps', 120); dt = duration/count
    state, field = initial.copy(), initial_field
    maximum_residual = 0.; history = []
    for step in range(count):
        penultimate, previous_field = state.copy(), field
        state, field, diagnostics = rk4(system, state, field, dt)
        maximum_residual = max(maximum_residual, diagnostics['max_residual'])
        if (step+1) % max(1, count//10) == 0 or step+1 == count:
            history.append(dict(step=step+1, time=(step+1)*dt, density_min=float(state[..., 0].min()),
                                energy_change=grid_change(d, system, ref, initial, initial_field, state, field)))
    change = grid_change(d, system, ref, initial, initial_field, state, field)
    result = dict(spec=spec, status='complete', duration=duration, dt=dt, steps=count,
                  initial=field_norms(spec, ref, d, initial, system, initial_field, 0., 12),
                  final=field_norms(spec, ref, d, state, system, field, duration, 12),
                  grid_change=change, max_residual=maximum_residual, history=history,
                  elapsed_seconds=time.monotonic()-start)
    if spec['kind'] == 'wave':
        result['energy_change_over_wave_energy'] = change['total']/ref.perturbation_energy()
    arrays = dict(initial_state=initial, final_state=state, penultimate_state=penultimate,
        initial_magnetic=initial_field.packed_increment(system.geometry),
        final_magnetic=field.packed_increment(system.geometry),
        penultimate_magnetic=previous_field.packed_increment(system.geometry), dt=np.array(dt))
    return result, arrays
