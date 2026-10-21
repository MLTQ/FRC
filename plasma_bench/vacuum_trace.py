"""Exterior magnetic traces for conservative interface stress and Poynting flux."""
import numpy as np
from scipy.constants import mu_0 as MU0


def tangential_field(mapping, total_flux, source=None):
    """Second-order vacuum-side derivatives at the actual gas interface.

    Surface-current force/work must use the exterior Maxwell traction; merely
    reflecting the interior tangential field silently loses interface work.
    With `source` (the analytic coil Field on the grid), `total_flux` must be the
    INDUCED flux only: its smooth harmonic exterior is differenced and the coil
    field is added analytically, so a coil filament inside the one-sided stencil
    cannot corrupt the trace.
    Returns (bz_wall, br_bottom, br_top, bz_inner); bz_inner is the inner-vacuum trace on the inner wall
    faces of an annular chamber (differenced one-sided toward the axis), else None.
    """
    m, mesh, grid = mapping, mapping.fluid, mapping.field
    i0 = m.r_offset
    i, lo, hi = i0+mesh.nr, m.z_offset, m.z_offset+mesh.nz
    if source is not None:
        bz_wall, br_bottom, br_top, bz_inner = tangential_field(mapping, total_flux)
        bz_nodes = source.bz[i, lo:hi+1]
        bz_wall = bz_wall+.5*(bz_nodes[:-1]+bz_nodes[1:])
        w = (mesh.radial_centroid-mesh.rf[:-1])/mesh.dr
        ends = [(1-w)*source.br[i0:i+1, j][:-1]+w*source.br[i0:i+1, j][1:] for j in (lo, hi)]
        if bz_inner is not None:
            inner_nodes = source.bz[i0, lo:hi+1]
            bz_inner = bz_inner+.5*(inner_nodes[:-1]+inner_nodes[1:])
        return bz_wall, br_bottom+ends[0], br_top+ends[1], bz_inner
    dr = (-3*total_flux[i]+4*total_flux[i+1]-total_flux[i+2])/(2*grid.dr)
    bz_nodes = dr[lo:hi+1]/mesh.radius
    bz_wall = .5*(bz_nodes[:-1]+bz_nodes[1:])
    bz_inner = None
    if i0:
        inner = (3*total_flux[i0]-4*total_flux[i0-1]+total_flux[i0-2])/(2*grid.dr)
        inner_nodes = inner[lo:hi+1]/mesh.inner_radius
        bz_inner = .5*(inner_nodes[:-1]+inner_nodes[1:])
    br_ends = []
    for j, direction in ((lo, -1), (hi, 1)):
        rows = slice(i0, i+1)
        dz = direction*(-3*total_flux[rows, j]+4*total_flux[rows, j+direction]-total_flux[rows, j+2*direction])/(2*grid.dz)
        br_nodes = np.divide(-dz, mesh.rf, out=np.zeros_like(dz), where=mesh.rf > 0)
        w = (mesh.radial_centroid-mesh.rf[:-1])/mesh.dr
        br_ends.append((1-w)*br_nodes[:-1]+w*br_nodes[1:])
    return bz_wall, br_ends[0], br_ends[1], bz_inner


def interface_resistivity(mapping, q, gamma, speed_cap=None):
    """Positive LLF-scale numerical resistance on plasma interface nodes only.

    A half-cell dual area gives E_num approximately a*jump(B_t)/2 for a
    surface current: D_num=a*h/4. The stiffness current supplies its sign.
    `speed_cap` (optional) bounds the local signal speed a; in near-vacuum wall
    cells a -> v_A grows without bound and so would this numerical dissipation.
    """
    mesh, grid = mapping.fluid, mapping.field
    if q.shape[-1] == 8:   # 2.5D: include B_theta = sqrt(<r^2>) g in the signal speed
        b2 = q[..., 4]**2+q[..., 5]**2+.5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]*q[..., 7]**2
        speed = np.hypot(q[..., 1], q[..., 2])+np.sqrt(gamma*q[..., 3]/q[..., 0]+b2/(MU0*q[..., 0]))
    else:
        speed = np.hypot(q[..., 1], q[..., 2])+np.sqrt(gamma*q[..., 3]/q[..., 0]+(q[..., 4]**2+q[..., 5]**2)/(MU0*q[..., 0]))
    if speed_cap is not None:
        speed = np.minimum(speed, speed_cap)
    mask = np.zeros(grid.shape)
    lo, hi = mapping.z_offset, mapping.z_offset+mesh.nz
    i0 = mapping.r_offset
    mask[i0+mesh.nr, lo:hi+1] = 1
    mask[max(i0, 1):i0+mesh.nr+1, [lo, hi]] = 1
    if i0:
        mask[i0, lo:hi+1] = 1
    return MU0*.25*min(mesh.dr, mesh.dz)*np.asarray(mapping.heat_to_cells.T@speed.ravel())*mapping.grid_vector(mask)
