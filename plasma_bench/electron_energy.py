"""Conservative pressure Poynting plus electron enthalpy energy transport."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .compatible_poynting import radial_faces
from .dual_enthalpy import enthalpy_flux
from .toroidal_state import radius_squared


def pressure_poynting_flux(model, g, radial_electric, axial_electric):
    """Integrated Epressure cross B power [W]; the poloidal cross term is toroidal."""
    geo, mesh = model.geometry, model.geometry.mesh
    er = np.asarray(radial_electric).reshape(mesh.nr, mesh.nz)
    ez = geo.scalar(axial_electric, 'pressure Ez')
    gf = radial_faces(mesh, g)
    gf[-1] = model.wall_btheta/mesh.radius
    radial = -2*np.pi*mesh.rf[:, None]**2*mesh.dz*ez*gf/MU0
    gz = .5*(g+np.roll(g, 1, axis=1))
    axial = 2*np.pi*geo.dx[:, None]*radius_squared(mesh)*er*gz/MU0
    return radial, np.concatenate((axial, axial[:, :1]), axis=1)


def energy_flux(model, g, entropy_density):
    _, _, er, ez = model.rhs(g, entropy_density)
    poynting = pressure_poynting_flux(model, g, er, ez)
    enthalpy = enthalpy_flux(model, g, entropy_density)
    return tuple(p+h for p, h in zip(poynting, enthalpy))


def energy_rate(model, g, entropy_density):
    fr, fz = energy_flux(model, g, entropy_density)
    return -(np.diff(fr, axis=0)+np.diff(fz, axis=1))/model.geometry.mesh.volumes
