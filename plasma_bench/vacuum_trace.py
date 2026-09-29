"""Exterior magnetic traces for conservative interface stress and Poynting flux."""
import numpy as np
from scipy.constants import mu_0 as MU0


def tangential_field(mapping, total_flux):
    """Second-order vacuum-side derivatives at the actual gas interface.

    Surface-current force/work must use the exterior Maxwell traction; merely
    reflecting the interior tangential field silently loses interface work.
    """
    m, mesh, grid = mapping, mapping.fluid, mapping.field
    i, lo, hi = mesh.nr, m.z_offset, m.z_offset+mesh.nz
    dr = (-3*total_flux[i]+4*total_flux[i+1]-total_flux[i+2])/(2*grid.dr)
    bz_nodes = dr[lo:hi+1]/mesh.radius
    bz_wall = .5*(bz_nodes[:-1]+bz_nodes[1:])
    br_ends = []
    for j, direction in ((lo, -1), (hi, 1)):
        dz = direction*(-3*total_flux[:i+1, j]+4*total_flux[:i+1, j+direction]-total_flux[:i+1, j+2*direction])/(2*grid.dz)
        br_nodes = np.divide(-dz, mesh.rf, out=np.zeros_like(dz), where=mesh.rf > 0)
        w = (mesh.radial_centroid-mesh.rf[:-1])/mesh.dr
        br_ends.append((1-w)*br_nodes[:-1]+w*br_nodes[1:])
    return bz_wall, br_ends[0], br_ends[1]


def interface_resistivity(mapping, q, gamma):
    """Positive LLF-scale numerical resistance on plasma interface nodes only.

    A half-cell dual area gives E_num approximately a*jump(B_t)/2 for a
    surface current: D_num=a*h/4. The stiffness current supplies its sign.
    """
    mesh, grid = mapping.fluid, mapping.field
    speed = np.hypot(q[..., 1], q[..., 2])+np.sqrt(gamma*q[..., 3]/q[..., 0]+(q[..., 4]**2+q[..., 5]**2)/(MU0*q[..., 0]))
    mask = np.zeros(grid.shape)
    lo, hi = mapping.z_offset, mapping.z_offset+mesh.nz
    mask[mesh.nr, lo:hi+1] = 1
    mask[1:mesh.nr+1, [lo, hi]] = 1
    return MU0*.25*min(mesh.dr, mesh.dz)*np.asarray(mapping.heat_to_cells.T@speed.ravel())*mapping.grid_vector(mask)
