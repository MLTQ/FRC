"""Independent physical reconstruction and parent-equation error norms."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from .current_transport_sampling import TransportSampler
from .regular_plasma_reference import exact


def norm(sampler, value, target, scale=None):
    error = value-target
    if error.ndim == sampler.weight.ndim+1:
        square = lambda x: np.sum(x*x, axis=-1)
        magnitude = lambda x: np.linalg.norm(x, axis=-1)
    else:
        square = lambda x: x*x
        magnitude = abs
    denominator = float(np.sum(sampler.weight*square(target))) if scale is None else scale**2*float(sampler.weight.sum())
    denominator = max(denominator, np.finfo(float).tiny)
    result = dict(relative_L2=float(np.sqrt(np.sum(sampler.weight*square(error))/denominator)),
                  max_absolute=float(np.max(magnitude(error))))
    for name, mask in [('axis_element', sampler.r <= sampler.d.radius/sampler.d.elements),
                       ('wall_element', sampler.r >= sampler.d.radius*(1-1/sampler.d.elements))]:
        local = float(np.sum(sampler.weight[mask]*square(target)[mask])) if scale is None else scale**2*float(sampler.weight[mask].sum())
        result[name] = float(np.sqrt(np.sum(sampler.weight[mask]*square(error)[mask])/max(local, np.finfo(float).tiny)))
    return result


def measure(spec, ref, d, u, system, field, rate, order=16):
    sample = TransportSampler(d, system.geometry.mesh, order)
    n, v, _, nt, vt = sample.fluid(u, rate.fluid)
    target = exact(ref, spec['kind'], sample.R, sample.Z)
    s = np.stack([sample.scalar(u[..., k])/n for k in (4, 5)], axis=-1)
    ct = np.stack([sample.scalar(rate.fluid[..., k]) for k in (4, 5)], axis=-1)
    h = E*np.exp(s)*(n[..., None]/d.thermal.n0)**(d.thermal.gammas-1)/(d.thermal.gammas-1)
    p = (d.thermal.gammas-1)*n[..., None]*h
    pt = p*(ct+(d.thermal.gammas-s)*nt[..., None])/n[..., None]
    b = sample.magnetic(field.packed_increment(system.geometry)); b[..., 2] += field.axial_bias
    bt = sample.magnetic(rate.magnetic)
    from plasma_bench.regular_current import RegularCurrent
    model = RegularCurrent(system.mapping, u[..., 0], ref.ion_mass)
    j = sample.electric(field.current(model)); jt = sample.electric(rate.current)
    fluid_power = d.mass*(.5*nt*np.sum(v*v, axis=-1)+n*np.sum(v*vt, axis=-1))
    fluid_power += np.sum(h*(ct+(d.thermal.gammas-s)*nt[..., None]), axis=-1)
    current_power = model.reduced_mass/E**2*(np.sum(j*jt, axis=-1)/n-.5*np.sum(j*j, axis=-1)*nt/n**2)
    magnetic_power = np.sum(b*bt, axis=-1)/MU0
    powers = {k: float(np.sum(sample.weight*a)) for k, a in
              [('fluid', fluid_power), ('current', current_power), ('magnetic', magnetic_power)]}
    velocity_scale = ref.velocity/ref.period if spec['kind'] == 'contact' else None
    norms = dict(density_rate=norm(sample, nt, target['nt']),
                 velocity_rate=norm(sample, vt, target['velocity_rate'], velocity_scale),
                 electron_pressure_rate=norm(sample, pt[..., 0], target['pressure_rate'][..., 0]),
                 ion_pressure_rate=norm(sample, pt[..., 1], target['pressure_rate'][..., 1]),
                 electric=norm(sample, sample.electric(rate.electric), target['electric']))
    if spec['kind'] == 'wave':
        norms.update(magnetic_rate=norm(sample, bt, target['magnetic_rate']),
                     current_rate=norm(sample, jt, target['current_rate']))
    else:
        norms['spurious_magnetic_rate_max'] = float(np.max(np.linalg.norm(bt, axis=-1)))
        norms['spurious_current_rate_max'] = float(np.max(np.linalg.norm(jt, axis=-1)))
    return dict(norms=norms, continuous_powers=powers, continuous_total_power=sum(powers.values()))
