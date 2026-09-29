"""Electron enthalpy flux through primal faces from a conservative dual current."""
import numpy as np
from scipy.constants import elementary_charge as E


def scalar_current_flux(geometry, radial, axial, scalar):
    """Integrate Q1 scalar times dual RT0 current over each primal face [scalar A]."""
    mesh = geometry.mesh
    value = geometry.scalar(scalar)
    radial, axial = np.asarray(radial, float), np.asarray(axial, float)
    if radial.shape != (mesh.nr+2, mesh.nz) or axial.shape != value.shape or not np.isfinite(radial).all() or not np.isfinite(axial).all():
        raise ValueError('Finite integrated dual current shapes required')
    edges = geometry.dual_x_faces
    width = np.diff(edges)
    fraction = ((geometry.x-edges[:-1])/width)[:, None]
    at_x = (1-fraction)*radial[:-1]+fraction*radial[1:]
    next_z = np.roll(value, -1, axis=1)
    fr = .5*(at_x*(.75*value+.25*next_z)+np.roll(at_x, -1, axis=1)*(.25*value+.75*next_z))
    at_z = .5*(axial+np.roll(axial, 1, axis=1))/width[:, None]
    fz = .5*geometry.dx[:, None]*(at_z[:-1]*(.75*value[:-1]+.25*value[1:])
                                 +at_z[1:]*(.25*value[:-1]+.75*value[1:]))
    return fr, np.concatenate((fz, fz[:, :1]), axis=1)


def enthalpy_flux(model, g, entropy_density):
    """Outward physical electron-current enthalpy flux, integrated watts."""
    temperature_j = model.pressure(entropy_density)/model.density
    current = model.geometry.current_flux(g, model.wall_btheta)
    flux = scalar_current_flux(model.geometry, *current, temperature_j)
    return tuple(-model.gamma/((model.gamma-1)*E)*f for f in flux)
