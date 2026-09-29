"""Positive face-field magnetic quadrature matching the variational edge energy."""
import numpy as np
from scipy.constants import mu_0 as MU0


def energy_density(mesh, psi):
    """Cell magnetic energy from squared staggered fields, before averaging."""
    br, bz = mesh.face_fields(psi)
    radial = (mesh.rf[:-1, None]*br[:-1]**2+mesh.rf[1:, None]*br[1:]**2)/(2*mesh.r[:, None])
    axial = .5*(bz[:, :-1]**2+bz[:, 1:]**2)
    return (radial+axial)/(2*MU0)


def selected_energy(mesh, psi, quadrature="cell"):
    if quadrature == "staggered":
        return energy_density(mesh, psi)
    if quadrature == "cell":
        br, bz = mesh.cell_fields(psi)
        return (br*br+bz*bz)/(2*MU0)
    raise ValueError("Magnetic energy quadrature must be cell or staggered")


def stiffness_action(mesh, psi):
    """Gradient of fluid-region edge energy, divided by 2 pi / mu0.

    Boundary edges have half dual measure on the fluid side. Axis psi is fixed;
    its returned reaction row is not an independent magnetic degree of freedom.
    """
    radial = mesh.dz*np.diff(psi, axis=0)/(mesh.r[:, None]*mesh.dr)
    radial[:, [0, -1]] *= .5
    axial = np.divide(mesh.dr*np.diff(psi, axis=1), mesh.rf[:, None]*mesh.dz,
                      out=np.zeros((mesh.nr+1, mesh.nz)), where=mesh.rf[:, None] > 0)
    axial[-1] *= .5
    result = np.zeros_like(psi)
    result[:-1] -= radial; result[1:] += radial
    result[:, :-1] -= axial; result[:, 1:] += axial
    return result
