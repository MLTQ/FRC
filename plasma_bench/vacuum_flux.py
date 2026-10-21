"""Conservative MHD fluxes with reflecting gas and transparent magnetic walls."""
import numpy as np
from scipy.constants import mu_0 as MU0
from .mhd_mesh import MHDMesh
from .mhd_reconstruction import Reconstruction, staggered_corner_fields
from .mhd_flux import rusanov, ghost_cells
from .mhd_induction import uct_electric
from .flux_update import FluxUpdate


def flux_components(mesh, q, psi, gamma, resistive_electric, boundary_field, interface_electric, spatial_order=2,
                    local_bounds=False, open_ends=False, magnetic_energy_flux="face", toroidal_electric=None,
                    toroidal_boundary=0.):
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
    An annular mesh (inner_radius > 0) has a second reflecting, magnetically transparent wall at the inner
    radius: boundary_field then carries a fourth entry, the inner-vacuum tangential Bz on the inner faces.

    2.5D (q with Omega = v_theta/r and g = B_theta/r columns): the gas also carries angular momentum and
    B_theta/r. Every gas boundary is insulating with an exact vacuum beyond it, so g is reflected odd there
    (B_theta = 0 on the face) and Omega even (slip). `toroidal_electric = (E_z, E_r, E_z_axis)` adds the
    non-ideal (resistive, Hall) E_z on radial faces and dr-averaged E_r on axial faces to the B_theta/r
    flux (r F_r = -E_z, F_z = E_r/r_c: Faraday's law on each r-z cell); E_z_axis enters through
    `FluxUpdate.toroidal_axis`. Their Poynting flux, and in the "dissipative"/"ct" modes that of the
    Rusanov B_theta/r flux, is paired with `toroidal_ohm.face_toroidal_field` (zero on every boundary without a
    threading current).
    `toroidal_boundary` (an annular chamber threaded by a current) is the vacuum r B_theta beyond the walls and
    ends: g is then reflected odd about its boundary value, and the boundary B_theta Poynting flux is returned
    separately in `FluxUpdate.toroidal_boundary` so it survives the variational boundary-flux substitution.
    """
    toroidal = q.shape[-1] == 8
    linked = toroidal and toroidal_boundary != 0
    r2_rows = .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]

    def odd(values, wall):
        """B_theta/r reflected odd about its value `wall` on the boundary face (zero without a linked current)."""
        values *= -1
        if linked:
            values += 2*wall
    if open_ends:
        # Outflow-only ends, per radial row: zero-gradient ghosts where the end cell flows
        # outward; the exact mirrored reflecting ghosts where it flows inward (no inflow).
        lo_wall, hi_wall = q[:, 1::-1].copy(), q[:, :-3:-1].copy()
        lo_wall[..., 2] *= -1; hi_wall[..., 2] *= -1
        lo_out, hi_out = np.repeat(q[:, :1], 2, axis=1), np.repeat(q[:, -1:], 2, axis=1)
        if toroidal:   # the vacuum's B_theta beyond either end, whether the gas flows out or not
            for ghost in (lo_wall, hi_wall, lo_out, hi_out):
                odd(ghost[..., 7], toroidal_boundary/r2_rows)
        lo = np.where((q[:, 0, 2] < 0)[:, None, None], lo_out, lo_wall)
        hi = np.where((q[:, -1, 2] > 0)[:, None, None], hi_out, hi_wall)
    else:
        lo, hi = q[:, 1::-1].copy(), q[:, :-3:-1].copy()
        lo[..., 2] *= -1; hi[..., 2] *= -1
        if toroidal:
            odd(lo[..., 7], toroidal_boundary/r2_rows); odd(hi[..., 7], toroidal_boundary/r2_rows)
    extended = np.concatenate((lo, q, hi), axis=1)
    ghost_mesh = MHDMesh(mesh.radius, mesh.half_length+2*mesh.dz, mesh.nr, mesh.nz+4, mesh.inner_radius)
    annular = mesh.inner_radius > 0
    reconstruction = (Reconstruction(ghost_mesh, extended, local_bounds=local_bounds,
                                     linked=toroidal_boundary if linked else 0., order=spatial_order)
                      if spatial_order == 2 or linked else None)
    if reconstruction is not None:
        left, right, down, up = reconstruction.faces()
        corners = list(reconstruction.corners())
        for index in (1, 3):
            corners[index][-1, :, 4] *= -1
            if toroidal:
                odd(corners[index][-1, :, 7], toroidal_boundary/mesh.rf[-1]**2)
        if annular:                     # the inner wall is a transparent wall, not the regular axis
            for index in (0, 2):
                corners[index][0, :, 4] *= -1
                if toroidal:
                    odd(corners[index][0, :, 7], toroidal_boundary/mesh.rf[0]**2)
        source = reconstruction.radial_source()
    else:
        ghost = ghost_cells(extended)
        ghost[-1, :, 4] *= -1
        if toroidal:
            odd(ghost[-1, :, 7], toroidal_boundary/mesh.rf[-1]**2)
        if annular:
            ghost[0, :, 4] *= -1
            if toroidal:
                odd(ghost[0, :, 7], toroidal_boundary/mesh.rf[0]**2)
        left, right = ghost[:-1, 1:-1].copy(), ghost[1:, 1:-1].copy()
        down, up = ghost[1:-1, :-1].copy(), ghost[1:-1, 1:].copy()
        corners = [ghost[:-1, :-1], ghost[1:, :-1], ghost[:-1, 1:], ghost[1:, 1:]]
        source = (extended[..., 3]+(extended[..., 4]**2+extended[..., 5]**2)/(2*MU0))/mesh.r[:, None]
        if toroidal:   # integral of (B_theta^2/2 mu0 + rho v_theta^2 - B_theta^2/mu0) dr / integral(r dr)
            source = source+(extended[..., 0]*extended[..., 6]**2-extended[..., 7]**2/(2*MU0))*mesh.radial_centroid[:, None]
    # The existing reconstruction reflects Br for a conducting radial wall.
    # Here the field is continuous through a magnetically transparent gas wall.
    if reconstruction is not None:
        right[-1, :, 4] *= -1
        if toroidal:
            odd(right[-1, :, 7], toroidal_boundary/mesh.rf[-1]**2)
        if annular:
            left[0, :, 4] *= -1
            if toroidal:
                odd(left[0, :, 7], toroidal_boundary/mesh.rf[0]**2)
    left, right, down, up = [value[:, 2:-2] for value in (left, right, down, up)]
    bz_wall, br_bottom, br_top = boundary_field[:3]
    left[-1, :, 5] = right[-1, :, 5] = bz_wall
    if annular:
        left[0, :, 5] = right[0, :, 5] = boundary_field[3]
    down[:, 0, 4] = up[:, 0, 4] = br_bottom
    down[:, -1, 4] = up[:, -1, 4] = br_top
    brf, bzf = mesh.face_fields(psi)
    mode = magnetic_energy_flux
    if toroidal:
        axial_radius = np.sqrt(.5*(mesh.rf[:-1]**2+mesh.rf[1:]**2))[:, None]
        fr, speed_r = rusanov(left, right, brf, 0, gamma, mode == "face", mode != "ct", radius=mesh.rf[:, None])
        fz, speed_z = rusanov(down, up, bzf, 1, gamma, mode == "face", mode != "ct", radius=axial_radius)
    else:
        fr, _ = rusanov(left, right, brf, 0, gamma, mode == "face", mode != "ct")
        fz, _ = rusanov(down, up, bzf, 1, gamma, mode == "face", mode != "ct")
    corners = [value[:, 2:-2] for value in corners]
    br_extended = np.pad(brf, ((0, 0), (2, 2)), mode="symmetric")
    bz_extended = np.pad(bzf, ((0, 0), (2, 2)), mode="reflect")
    if spatial_order == 2:
        fields = staggered_corner_fields(ghost_mesh, br_extended, bz_extended)
    else:
        bz_pad = np.pad(bz_extended, ((1, 1), (0, 0)), mode="edge")
        br_pad = np.concatenate((br_extended[:, -1:], br_extended, br_extended[:, :1]), axis=1)
        fields = (bz_pad[:-1], bz_pad[1:], br_pad[:, :-1], br_pad[:, 1:])
    fields = [value[:, 2:-2] for value in fields]
    central, dissipative = uct_electric(corners, fields, gamma, parts=True,
                                        radius=mesh.rf[:, None] if toroidal else None)
    electric = central+dissipative+resistive_electric+interface_electric
    if not annular:
        electric[0] = 0  # Regular axis; no artificial zero at the chamber walls/ends.
        dissipative[0] = 0
        central[0] = 0
    nonideal = resistive_electric+interface_electric
    if mode in ("dissipative", "ct"):
        nonideal = nonideal+dissipative
    if mode == "ct":
        nonideal = nonideal+central
    fr[..., 3] += .5*(nonideal[:, :-1]+nonideal[:, 1:])*.5*(left[..., 5]+right[..., 5])/MU0
    fz[..., 3] -= .5*(nonideal[:-1]+nonideal[1:])*.5*(down[..., 4]+up[..., 4])/MU0
    axis_electric = boundary_poynting = None
    if toroidal:
        axis_electric, boundary_poynting = _toroidal(mesh, q, fr, fz, speed_r, speed_z, left, right, down, up,
                                                     toroidal_electric, mode, toroidal_boundary)
    result = FluxUpdate(fr, fz, electric, source[:, 2:-2], toroidal_axis=axis_electric,
                        toroidal_boundary=boundary_poynting)
    result.dissipative_electric = dissipative if mode != "face" else None
    return result


def _toroidal(mesh, q, fr, fz, speed_r, speed_z, left, right, down, up, toroidal_electric, mode, boundary=0.):
    """Add the non-ideal B_theta/r flux and the toroidal Poynting flux; return the axis E_z (or None) and, with
    a linked vacuum r B_theta (`boundary`), the boundary Poynting fluxes (else None)."""
    from .toroidal_ohm import face_toroidal_field
    ez, er, ez_axis = toroidal_electric
    rf, rc = mesh.rf[:, None], mesh.r[:, None]
    # Electric fields that the Rusanov B_theta/r flux represents: E_z = -r F_r, dr-averaged E_r = r_c F_z.
    if mode == "ct":
        ez_poynting, er_poynting = -rf*fr[..., 5]+ez, rc*fz[..., 5]+er
    elif mode == "dissipative":
        ez_poynting = .5*rf*speed_r*(right[..., 7]-left[..., 7])+ez
        er_poynting = -.5*rc*speed_z*(up[..., 7]-down[..., 7])+er
    else:
        ez_poynting, er_poynting = ez, er
    first = 0 if mesh.annular else 1    # the axis face (r = 0) takes E_z through FluxUpdate.toroidal_axis
    fr[first:, :, 5] -= ez[first:]/rf[first:]
    fz[..., 5] += er/rc
    b_radial, b_axial = face_toroidal_field(mesh, .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]*q[..., 7], boundary)
    fr[..., 3] -= ez_poynting*b_radial/MU0
    fz[..., 3] += er_poynting*b_axial/MU0
    edges = None
    if boundary:
        edges = (-ez_poynting[-1]*b_radial[-1]/MU0, er_poynting[:, 0]*b_axial[:, 0]/MU0,
                 er_poynting[:, -1]*b_axial[:, -1]/MU0, -ez_poynting[0]*b_radial[0]/MU0 if mesh.annular else None)
    return None if mesh.annular else ez_axis, edges


def flux_rhs(mesh, q, psi, gamma, resistive_electric, boundary_field, interface_electric, spatial_order=2):
    return flux_components(mesh, q, psi, gamma, resistive_electric, boundary_field,
                           interface_electric, spatial_order).rhs(mesh)[:3]
