"""Independent physical values and stable continuous conservation changes."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from plasma_bench.regular_current import RegularCurrent


def values(sample, state, magnetic, system, bias):
    n, velocity, _ = sample.fluid(state)
    entropy = np.stack([sample.scalar(state[..., k])/n for k in (4, 5)], axis=-1)
    h = E*np.exp(entropy)*(n[..., None]/sample.d.thermal.n0)**(sample.d.thermal.gammas-1)/(sample.d.thermal.gammas-1)
    model = RegularCurrent(system.mapping, state[..., 0], system.ion_mass)
    current = sample.electric(model.current(magnetic))
    field = sample.magnetic(magnetic)
    return dict(n=n, velocity=velocity, entropy=entropy, h=h, current=current,
                magnetic=field, reduced_mass=model.reduced_mass)


def changes(sample, system, bias, initial_state, initial_magnetic, state, magnetic):
    a = values(sample, initial_state, initial_magnetic, system, bias)
    b = values(sample, state, magnetic, system, bias)
    n0, n = a['n'], b['n']; difference = state-initial_state
    dn = sample.scalar(difference[..., 0])
    dc = np.stack([sample.scalar(difference[..., k]) for k in (4, 5)], axis=-1)
    ds = (dc-a['entropy']*dn[..., None])/n[..., None]
    thermal = n0[..., None]*a['h']*np.expm1(ds+sample.d.thermal.gammas*np.log1p(dn/n0)[..., None])
    kinetic = .5*sample.d.mass*(n*np.sum(b['velocity']**2, axis=-1)-n0*np.sum(a['velocity']**2, axis=-1))
    current = .5*b['reduced_mass']/E**2*(np.sum(b['current']**2, axis=-1)/n-np.sum(a['current']**2, axis=-1)/n0)
    field_difference = sample.magnetic(magnetic-initial_magnetic)
    field = np.sum(field_difference*(b['magnetic']+a['magnetic']), axis=-1)/(2*MU0)
    cross = bias*field_difference[..., 2]/MU0
    terms = dict(bulk_kinetic=kinetic, electron_internal=thermal[..., 0],
        ion_internal=thermal[..., 1], current=current, magnetic_increment=field, magnetic_bias_cross=cross)
    integrals = {name: float(np.sum(sample.weight*value)) for name, value in terms.items()}
    full_b = a['magnetic'].copy(); full_b[..., 2] += bias
    initial_energy = float(np.sum(sample.weight*(.5*sample.d.mass*n0*np.sum(a['velocity']**2, axis=-1)
        +n0*a['h'].sum(axis=-1)+.5*a['reduced_mass']/E**2*np.sum(a['current']**2, axis=-1)/n0
        +np.sum(full_b**2, axis=-1)/(2*MU0))))
    particles = float(np.sum(sample.weight*n0)); temperature = b['h']*(sample.d.thermal.gammas-1)/E
    return dict(energy_changes=integrals, energy_change=sum(integrals.values()),
        initial_energy=initial_energy, relative_energy_drift=sum(integrals.values())/initial_energy,
        particle_relative_drift=float(np.sum(sample.weight*dn)/particles),
        entropy_change_per_particle=[float(np.sum(sample.weight*dc[..., k])/particles) for k in range(2)],
        point_density_min=float(state[..., 0].min()), sampled_density_min=float(n.min()),
        temperature_min=temperature.min(axis=(0, 1)).tolist(), temperature_max=temperature.max(axis=(0, 1)).tolist())


def signal(sample, system, negative_state, positive_state, negative_magnetic, positive_magnetic, epsilon):
    minus = values(sample, negative_state, negative_magnetic, system, 0.)
    plus = values(sample, positive_state, positive_magnetic, system, 0.)
    delta = positive_state-negative_state; dn = sample.scalar(delta[..., 0])
    dc = np.stack([sample.scalar(delta[..., k]) for k in (4, 5)], axis=-1)
    ds = (dc-minus['entropy']*dn[..., None])/plus['n'][..., None]
    physical = np.concatenate((dn[..., None], plus['velocity']-minus['velocity'], ds,
        sample.magnetic(positive_magnetic-negative_magnetic)), axis=-1)/(2*epsilon)
    pressure = minus['n'][..., None]*minus['h']*(sample.d.thermal.gammas-1)
    pressure *= np.expm1(ds+sample.d.thermal.gammas*np.log1p(dn/minus['n'])[..., None])/(2*epsilon)
    return physical, pressure, (plus['current']-minus['current'])/(2*epsilon)
