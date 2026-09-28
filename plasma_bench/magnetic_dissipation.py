"""Explicitly measured, vanishing magnetic diffusion for compressive shocks."""
import numpy as np
from scipy.constants import mu_0 as MU0


def shock_resistance(mapping, primitive, total_flux, velocity, coefficient, gamma):
    """Return numerical R and eta; eta_num=mu0 C h min(h compression, speed).

    This is numerical stabilization, not a material transport coefficient.
    It vanishes as h² in smooth compressive flow and as h across a shock.
    """
    if coefficient == 0:
        return np.zeros(mapping.node_count), 0.
    mesh = mapping.fluid
    vr, vz = [v.reshape(mesh.nr, mesh.nz) for v in np.split(velocity, 2)]
    divergence = np.gradient(vr, mesh.radial_centroid, axis=0, edge_order=2)+vr/mesh.radial_centroid[:, None]
    divergence += np.gradient(vz, mesh.z, axis=1, edge_order=2)
    br, bz = mesh.cell_fields(mapping.fluid_flux(total_flux))
    speed = np.hypot(vr, vz)+np.sqrt(gamma*primitive[..., 3]/primitive[..., 0]+(br*br+bz*bz)/(MU0*primitive[..., 0]))
    h = min(mesh.dr, mesh.dz)
    diffusivity = coefficient*h*np.minimum(h*np.maximum(-divergence, 0.), speed)
    # The transpose heat map averages adjoining fluid cells with nonnegative
    # weights; no resistivity is assigned to vacuum or multipole unknowns.
    eta = MU0*np.asarray(mapping.heat_to_cells.T@diffusivity.ravel())
    resistance = np.divide(eta*mapping.node_radius, MU0*mapping.node_area,
                            out=np.zeros(mapping.node_count), where=mapping.active)
    return resistance, float(eta.max(initial=0.))
