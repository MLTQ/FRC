"""Conservative MHD fluxes with reflecting gas and transparent magnetic walls."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd_mesh import MHDMesh
from .mhd_reconstruction import Reconstruction, staggered_corner_fields
from .mhd_flux import rusanov, ghost_cells
from .mhd_induction import uct_electric
from .flux_update import FluxUpdate


def flux_components(mesh, q, psi, gamma, resistive_electric, boundary_field, interface_electric, spatial_order=2,
                    local_bounds=False, open_ends=False, magnetic_energy_flux="face"):
    """Shared reconstructed states for fluid fluxes and constrained induction.

    The reflected gas traces are boundary data, not exterior fluid unknowns.
    Normal velocity reflects; magnetic normal field remains the shared CT field.
    magnetic_energy_flux: "face" (historical) transports magnetic energy with the face states:
    central (B^2/mu0)v_n - (v.B)B_n/mu0 plus the LLF jump of B_t^2/2mu0. "dissipative" replaces
    that LLF jump by the Poynting flux of the UCT dissipative electric field (the one that changes
    psi), paired exactly like resistivity. "ct" carries ALL magnetic energy transport as the
    Poynting flux E_theta x B_t/mu0 of the full CT electric field (ideal, dissipative, resistive,
    interface), so a current-free cell exchanges no spurious energy with its gas. The returned
    FluxUpdate has `dissipative_electric` set unless the mode is "face".
    """
    if open_ends:
        # Outflow-only ends, per radial row: zero-gradient ghosts where the end cell flows
        # outward; the exact mirrored reflecting ghosts where it flows inward (no inflow).
        lo_wall, hi_wall = q[:, 1::-1].copy(), q[:, :-3:-1].copy()
        lo_wall[..., 2] *= -1; hi_wall[..., 2] *= -1
        lo_out, hi_out = np.repeat(q[:, :1], 2, axis=1), np.repeat(q[:, -1:], 2, axis=1)
        lo = np.where((q[:, 0, 2] < 0)[:, None, None], lo_out, lo_wall)
        hi = np.where((q[:, -1, 2] > 0)[:, None, None], hi_out, hi_wall)
    else:
        lo, hi = q[:, 1::-1].copy(), q[:, :-3:-1].copy()
        lo[..., 2] *= -1; hi[..., 2] *= -1
    extended = np.concatenate((lo, q, hi), axis=1)
    ghost_mesh = MHDMesh(mesh.radius, mesh.half_length+2*mesh.dz, mesh.nr, mesh.nz+4)
    reconstruction = Reconstruction(ghost_mesh, extended, local_bounds=local_bounds) if spatial_order == 2 else None
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
    mode = magnetic_energy_flux
    fr, _ = rusanov(left, right, brf, 0, gamma, mode == "face", mode != "ct")
    fz, _ = rusanov(down, up, bzf, 1, gamma, mode == "face", mode != "ct")
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
    central, dissipative = uct_electric(corners, fields, gamma, parts=True)
    electric = central+dissipative+resistive_electric+interface_electric
    electric[0] = 0  # Regular axis; no artificial zero at the chamber wall/ends.
    dissipative[0] = 0
    central[0] = 0
    nonideal = resistive_electric+interface_electric
    if mode in ("dissipative", "ct"):
        nonideal = nonideal+dissipative
    if mode == "ct":
        nonideal = nonideal+central
    fr[..., 3] += .5*(nonideal[:, :-1]+nonideal[:, 1:])*.5*(left[..., 5]+right[..., 5])/MU0
    fz[..., 3] -= .5*(nonideal[:-1]+nonideal[1:])*.5*(down[..., 4]+up[..., 4])/MU0
    result = FluxUpdate(fr, fz, electric, source[:, 2:-2])
    result.dissipative_electric = dissipative if mode != "face" else None
    return result


def flux_rhs(mesh, q, psi, gamma, resistive_electric, boundary_field, interface_electric, spatial_order=2):
    return flux_components(mesh, q, psi, gamma, resistive_electric, boundary_field,
                           interface_electric, spatial_order).rhs(mesh)
