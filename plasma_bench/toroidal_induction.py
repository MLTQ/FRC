"""Toroidal-aware UCT and dissipative toroidal curl with its energy pairing."""
from functools import lru_cache
import numpy as np
from scipy.linalg import solve_banded
from scipy.constants import mu_0 as MU0
from .mhd_reconstruction import staggered_corner_fields
from .toroidal_flux import fast_speed
from .toroidal_state import radius_squared


def ideal_emf(mesh, reconstruction, br, bz, gamma, order):
    corners = reconstruction.corners()
    sw, se, nw, ne = corners
    if order == 2:
        bl, brt, bd, bu = staggered_corner_fields(mesh, br, bz)
    else:
        pad = np.pad(bz, ((1, 1), (0, 0)), mode='edge')
        bl, brt = pad[:-1], pad[1:]
        pad = np.concatenate((br[:, -1:], br, br[:, :1]), axis=1)
        bd, bu = pad[:, :-1], pad[:, 1:]
    r = mesh.rf[:, None]
    ar = np.maximum.reduce([abs(q[..., 1])+fast_speed(q, r, 0, gamma) for q in corners])
    az = np.maximum.reduce([abs(q[..., 2])+fast_speed(q, r, 1, gamma) for q in corners])
    electric = .25*((sw[..., 1]+nw[..., 1])*bl+(se[..., 1]+ne[..., 1])*brt
                    -(sw[..., 2]+se[..., 2])*bd-(nw[..., 2]+ne[..., 2])*bu)
    electric -= .5*ar*(brt-bl)
    electric += .5*az*(bu-bd)
    electric[[0, -1]] = 0.
    electric[:, -1] = electric[:, 0]
    return electric



@lru_cache(maxsize=16)
def radial_current_mass(nr, dr):
    """Positive compatible current mass, including regular-axis moments."""
    edge = (np.arange(nr)+.5)*dr**2/8  # Includes the prescribed-zero wall edge.
    weight = np.arange(nr)*dr**2
    weight[0] = dr**2/4
    adjacent = edge.copy()
    adjacent[1:] += edge[:-1]
    diagonal = weight-adjacent
    band = np.zeros((3, nr))
    band[0, 1:] = edge[:-1]
    band[1] = diagonal
    band[2, :-1] = edge[:-1]
    node_weight = weight-2*adjacent
    for array in (band, edge, node_weight):
        array.setflags(write=False)
    return band, edge, node_weight


def toroidal_resistivity(mesh, g, eta):
    """Return dg/dt, radial/axial Poynting flux, positive physical heat rate."""
    ez = np.zeros((mesh.nr+1, mesh.nz))
    er_over_r = np.zeros((mesh.nr, mesh.nz+1))
    if eta == 0:
        return np.zeros_like(g), ez, er_over_r, 0.
    d = eta/MU0
    s = radius_squared(mesh)
    band, edge, node_weight = radial_current_mass(mesh.nr, mesh.dr)
    weak_curl = s*g
    weak_curl = np.concatenate((weak_curl[:1], np.diff(weak_curl, axis=0)), axis=0)
    ez[:-1] = d*solve_banded((1, 1), band, weak_curl)
    # Ez(R)=0 is the conducting wall. Ez(axis) is solved, not set to zero.
    er_over_r[:, :-1] = -d*(g-np.roll(g, 1, axis=1))/mesh.dz
    er_over_r[:, -1] = er_over_r[:, 0]
    derivative = np.diff(ez, axis=0)/mesh.radial_measure[:, None]-np.diff(er_over_r, axis=1)/mesh.dz
    pad = np.pad(g, ((1, 1), (0, 0)), mode='edge')
    gr = .5*(pad[:-1]+pad[1:])
    gz = .5*(np.concatenate((g[:, -1:], g), axis=1)+np.concatenate((g, g[:, :1]), axis=1))
    sr = -ez*mesh.rf[:, None]*gr/MU0
    sz = er_over_r*s*gz/MU0
    radial_heat = 2*np.pi*mesh.dz/eta*(np.sum(node_weight[:, None]*ez[:-1]**2)
                  +np.sum(edge[:, None]*(ez[:-1]+ez[1:])**2))
    axial_heat = np.sum(mesh.volumes*s*er_over_r[:, :-1]**2)/eta
    return derivative, sr, sz, float(radial_heat+axial_heat)
