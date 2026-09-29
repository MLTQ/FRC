"""Conservative MHD fluxes with reflecting gas and transparent magnetic walls."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd_mesh import MHDMesh
from .mhd_reconstruction import Reconstruction, staggered_corner_fields
from .mhd_flux import rusanov, ghost_cells
from .mhd_induction import uct_electric
from .flux_update import FluxUpdate


def flux_components(mesh, q, psi, gamma, resistive_electric, boundary_field, interface_electric, spatial_order=2):
    """Shared reconstructed states for fluid fluxes and constrained induction.

    The reflected gas traces are boundary data, not exterior fluid unknowns.
    Normal velocity reflects; magnetic normal field remains the shared CT field.
    """
    lo, hi = q[:, 1::-1].copy(), q[:, :-3:-1].copy()
    lo[..., 2] *= -1; hi[..., 2] *= -1
    extended = np.concatenate((lo, q, hi), axis=1)
    ghost_mesh = MHDMesh(mesh.radius, mesh.half_length+2*mesh.dz, mesh.nr, mesh.nz+4)
    reconstruction = Reconstruction(ghost_mesh, extended) if spatial_order == 2 else None
    if reconstruction is not None:
        left, right, down, up = reconstruction.faces()
        corners = list(reconstruction.corners())
        for index in (1, 3):
            corners[index][-1, :, 4] *= -1
        source = reconstruction.radial_source()
    else:
        ghost = ghost_cells(extended)
        ghost[-1, :, 4] *= -1
        left, right = ghost[:-1, 1:-1].copy(), ghost[1:, 1:-1].copy()
        down, up = ghost[1:-1, :-1].copy(), ghost[1:-1, 1:].copy()
        corners = [ghost[:-1, :-1], ghost[1:, :-1], ghost[:-1, 1:], ghost[1:, 1:]]
        source = (extended[..., 3]+(extended[..., 4]**2+extended[..., 5]**2)/(2*MU0))/mesh.r[:, None]
    # The existing reconstruction reflects Br for a conducting radial wall.
    # Here the field is continuous through a magnetically transparent gas wall.
    if reconstruction is not None:
        right[-1, :, 4] *= -1
    left, right, down, up = [value[:, 2:-2] for value in (left, right, down, up)]
    bz_wall, br_bottom, br_top = boundary_field
    left[-1, :, 5] = right[-1, :, 5] = bz_wall
    down[:, 0, 4] = up[:, 0, 4] = br_bottom
    down[:, -1, 4] = up[:, -1, 4] = br_top
    brf, bzf = mesh.face_fields(psi)
    fr, _ = rusanov(left, right, brf, 0, gamma)
    fz, _ = rusanov(down, up, bzf, 1, gamma)
    corners = [value[:, 2:-2] for value in corners]
    br_extended = np.pad(brf, ((0, 0), (2, 2)), mode="symmetric")
    bz_extended = np.pad(bzf, ((0, 0), (2, 2)), mode="reflect")
    if reconstruction is not None:
        fields = staggered_corner_fields(ghost_mesh, br_extended, bz_extended)
    else:
        bz_pad = np.pad(bz_extended, ((1, 1), (0, 0)), mode="edge")
        br_pad = np.concatenate((br_extended[:, -1:], br_extended, br_extended[:, :1]), axis=1)
        fields = (bz_pad[:-1], bz_pad[1:], br_pad[:, :-1], br_pad[:, 1:])
    fields = [value[:, 2:-2] for value in fields]
    electric = uct_electric(corners, fields, gamma)+resistive_electric+interface_electric
    electric[0] = 0  # Regular axis; no artificial zero at the chamber wall/ends.
    nonideal = resistive_electric+interface_electric
    fr[..., 3] += .5*(nonideal[:, :-1]+nonideal[:, 1:])*.5*(left[..., 5]+right[..., 5])/MU0
    fz[..., 3] -= .5*(nonideal[:-1]+nonideal[1:])*.5*(down[..., 4]+up[..., 4])/MU0
    return FluxUpdate(fr, fz, electric, source[:, 2:-2])


def flux_rhs(mesh, q, psi, gamma, resistive_electric, boundary_field, interface_electric, spatial_order=2):
    return flux_components(mesh, q, psi, gamma, resistive_electric, boundary_field,
                           interface_electric, spatial_order).rhs(mesh)
