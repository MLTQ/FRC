"""Advection of electron entropy proxy with the actual dual bulk particle flux."""
import numpy as np
from .mhd_reconstruction import radial_slopes, axial_increments
from .mhd_positivity import positive_scale


def entropy_rate(mapping, density, entropy, radial, axial, *, order=2):
    """Return -div(Gamma * upwind(U/n)); no physical heating or drift term."""
    geo = mapping.geometry
    n, u = geo.scalar(density, 'number density'), geo.scalar(entropy, 'entropy proxy')
    radial, axial = np.asarray(radial, float), np.asarray(axial, float)
    if radial.shape != (mapping.mesh.nr+2, mapping.mesh.nz) or axial.shape != u.shape:
        raise ValueError('Dual integrated particle flux shapes required')
    if not np.isfinite(radial).all() or not np.isfinite(axial).all() or np.any(radial[[0, -1]] != 0):
        raise ValueError('Finite dual fluxes with impermeable radial boundaries required')
    if np.any(n <= 0) or np.any(u <= 0) or isinstance(order, bool) or order not in (1, 2):
        raise ValueError('Positive density/entropy and order 1 or 2 required')
    w = u/n
    low, high, down, up = (w.copy() for _ in range(4))
    if order == 2:
        slope = radial_slopes(w, mapping.centers_x)
        increment = axial_increments(w)
        lo = slope*(geo.dual_x_faces[:-1, None]-mapping.centers_x)
        hi = slope*(geo.dual_x_faces[1:, None]-mapping.centers_x)
        factor = positive_scale(w, w+np.minimum.reduce((lo, hi, -.5*abs(increment))))
        low, high = w+factor*lo, w+factor*hi
        down, up = w-.5*factor*increment, w+.5*factor*increment
    flux_r = np.zeros_like(radial)
    flux_r[1:-1] = radial[1:-1]*np.where(radial[1:-1] >= 0, high[:-1], low[1:])
    flux_z = axial*np.where(axial >= 0, up, np.roll(down, -1, axis=1))
    return -geo.divergence(flux_r, flux_z)
