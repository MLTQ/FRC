"""Full electric potential and conservative current/enthalpy/Poynting face power."""
import numpy as np
from scipy.constants import elementary_charge as E
from .compatible_poynting import energy_flux
from .dual_enthalpy import enthalpy_flux, scalar_current_flux


def current_particle_energy(model, current):
    """Nodal mu |w|²/2 [J]; declared traces for the kinetic flux/potential."""
    mesh, geo = model.mesh, model.electron.geometry
    jr, jt, jz = model.geometry.split_electric(current)
    jr = jr.reshape(mesh.nr, mesh.nz); jz = jz.reshape(mesh.nr+1, mesh.nz)
    radial = np.empty_like(jz); radial[0] = 0.; radial[-1] = 0.
    dx = geo.dx
    radial[1:-1] = mesh.rf[1:-1, None]*(dx[:-1, None]*jr[:-1]+dx[1:, None]*jr[1:])/(dx[:-1, None]+dx[1:, None])
    theta = np.pad(jt.reshape(mesh.nr-1, mesh.nz), ((1, 1), (0, 0)))
    axial = .5*(jz+np.roll(jz, 1, axis=1))
    return .5*model.reduced_mass*(radial*radial+theta*theta+axial*axial)/(E*model.electron.density)**2


def flux(model, packed, current, ue, ui, electric_reduced):
    """Return full E coefficients and integrated outward face powers [W]."""
    psi, g = model.geometry.unpack(packed)
    _, pressure_potential = model.pressure(ue, ui)
    kinetic_particle = current_particle_energy(model, current)
    electric = electric_reduced+model.gradient(pressure_potential-model.delta*kinetic_particle/E)
    sr, sz = energy_flux(model.geometry, psi, g, electric)
    mesh = model.mesh
    power = [2*np.pi*mesh.rf[:, None]*mesh.dz*sr, 2*np.pi*mesh.radial_measure[:, None]*sz]
    electron = enthalpy_flux(model.electron, g, ue); ion = enthalpy_flux(model.ion, g, ui)
    charge = model.electron.geometry.current_flux(g, 0.)
    drift = scalar_current_flux(model.electron.geometry, *charge, kinetic_particle)
    for k in range(2):
        power[k] += model.beta*electron[k]-model.alpha*ion[k]-model.delta/E*drift[k]
    return electric, tuple(power)
