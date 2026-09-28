"""Conservative gas-only Euler step for the split plasma/vacuum solver."""
import numpy as np
from .mhd_mesh import MHDMesh
from .mhd_state import MHDState, primitive, from_primitive
from .mhd_reconstruction import Reconstruction
from .mhd_flux import rusanov


def gas_primitive(mesh, conserved, gamma):
    return primitive(mesh, MHDState(conserved, np.zeros((mesh.nr+1, mesh.nz+1))), gamma)


def gas_from_primitive(mesh, density, pressure, *, vr=0., vz=0., gamma=5/3):
    return from_primitive(mesh, np.zeros((mesh.nr+1, mesh.nz+1)), density, pressure, vr, vz, gamma).conserved


def _rhs(mesh, conserved, gamma):
    q = gas_primitive(mesh, conserved, gamma)
    # Two reflected layers let the existing reconstruction operate on physical
    # faces without any periodic connection between the actual chamber ends.
    left, right = q[:, 1::-1].copy(), q[:, :-3:-1].copy()
    left[..., 2] *= -1; right[..., 2] *= -1
    extended = np.concatenate((left, q, right), axis=1)
    ghost_mesh = MHDMesh(mesh.radius, mesh.half_length+2*mesh.dz, mesh.nr, mesh.nz+4)
    reconstruction = Reconstruction(ghost_mesh, extended)
    radial_left, radial_right, axial_down, axial_up = reconstruction.faces()
    fr, _ = rusanov(radial_left[:, 2:-2], radial_right[:, 2:-2], 0., 0, gamma)
    fz, _ = rusanov(axial_down[:, 2:-2], axial_up[:, 2:-2], 0., 1, gamma)
    result = -np.diff(mesh.rf[:, None, None]*fr, axis=0)/mesh.radial_measure[:, None, None]
    result -= np.diff(fz, axis=1)/mesh.dz
    result[..., 1] += reconstruction.radial_source()[:, 2:-2]
    return result


def hydro_step(mesh, conserved, dt, gamma):
    """SSP RK2 with reflecting slip walls; caller supplies a valid CFL step."""
    a = _rhs(mesh, conserved, gamma)
    predictor = conserved+dt*a
    gas_primitive(mesh, predictor, gamma)
    b = _rhs(mesh, predictor, gamma)
    final = conserved+.5*dt*(a+b)
    gas_primitive(mesh, final, gamma)
    return final
