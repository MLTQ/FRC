"""Optional compressive bulk viscosity with matched mechanical work and heat."""
import numpy as np


def bulk_viscosity(fluid, state, fields, *, linear=.5, quadratic=2.):
    """Return a numerical bulk pressure, its nodal force, and positive heating."""
    if not np.isfinite([linear, quadratic]).all() or min(linear, quadratic) < 0:
        raise ValueError('Finite nonnegative numerical viscosity coefficients required')
    velocity = fields['velocity'][fluid.ids]
    divergence = np.einsum('cqna,cna->cq', fields['divergence'], velocity)
    compression = np.maximum(-divergence, 0.)
    # Smallest physical cell altitude, using all corner areas and edge lengths.
    x = state.coordinates[fluid.ids]; edge = np.roll(x, -1, axis=1)-x
    nxt = np.roll(edge, -1, axis=1)
    area = edge[..., 0]*nxt[..., 1]-edge[..., 1]*nxt[..., 0]
    altitude = np.min(area/np.maximum(np.linalg.norm(edge, axis=-1), np.linalg.norm(nxt, axis=-1)), axis=1)
    if np.any(altitude <= 0):
        raise ValueError('Positive cell altitude required; no viscous mesh repair')
    density = fluid.mass*fields['n']
    sound = np.sqrt(np.sum(fluid.gammas*fields['pressure'], axis=-1)/density)
    speed = altitude[:, None]*compression
    pressure = density*(linear*sound*speed+quadratic*speed**2)
    heating = pressure*compression
    local = np.einsum('cq,cqna->cna', fields['volume']*pressure, fields['divergence'])
    force = np.zeros_like(state.coordinates)
    np.add.at(force, fluid.ids.ravel(), local.reshape(-1, 2)); force[fluid.domain.axis, 0] = 0.
    power = float(np.sum(fields['volume']*heating))
    mechanical = float(np.sum(force*fields['velocity']))
    return dict(force=force, heating=heating, pressure=pressure, divergence=divergence,
                power=power, mechanical_power=mechanical, power_residual=mechanical+power,
                minimum_cell_altitude_m=float(altitude.min()))
