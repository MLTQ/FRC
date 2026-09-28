"""Compatible flux reconstruction and the axisymmetric diffusion stiffness."""
import numpy as np
from scipy.sparse import coo_matrix
from .geometry import Grid


def field_from_flux(grid: Grid, psi):
    psi = np.asarray(psi)
    if psi.shape != grid.shape:
        raise ValueError("Flux must have grid shape")
    r = grid.r[:, None]
    pr, pz = np.gradient(psi, grid.dr, grid.dz, edge_order=2)
    br = np.divide(-pz, r, out=np.zeros_like(pz), where=r > 0)
    bz = np.divide(pr, r, out=np.zeros_like(pr), where=r > 0)
    bz[0] = (8*psi[1]-0.5*psi[2])/(3*grid.dr**2)
    return br, bz


def cell_divergence(grid: Grid, psi):
    """Divergence of face fluxes; exactly telescopes, including axis cells."""
    radial_flux = -np.diff(psi, axis=1)  # r Br dz
    axial_flux = np.diff(psi, axis=0)    # integral r Bz dr
    net = np.diff(radial_flux, axis=0)+np.diff(axial_flux, axis=1)
    volumes = (np.diff(grid.r**2)/2)[:, None]*grid.dz
    return net/volumes


def stiffness(grid: Grid):
    """SPD interior-node energy matrix; homogeneous flux on all boundaries."""
    ni, nj = grid.nr-2, grid.nz-2
    rows, cols, data = [], [], []
    for i in range(ni):
        r = grid.r[i+1]
        for j in range(nj):
            p = i*nj+j
            diagonal = 0.0
            for di, dj, weight in ((-1, 0, grid.dz/((r-grid.dr/2)*grid.dr)),
                                    (1, 0, grid.dz/((r+grid.dr/2)*grid.dr)),
                                    (0, -1, grid.dr/(r*grid.dz)),
                                    (0, 1, grid.dr/(r*grid.dz))):
                diagonal += weight
                if 0 <= i+di < ni and 0 <= j+dj < nj:
                    rows.append(p); cols.append((i+di)*nj+j+dj); data.append(-weight)
            rows.append(p); cols.append(p); data.append(diagonal)
    return coo_matrix((data, (rows, cols)), shape=(ni*nj, ni*nj)).tocsc()


def expand(grid: Grid, interior):
    out = np.zeros(grid.shape)
    out[1:-1, 1:-1] = np.reshape(interior, (grid.nr-2, grid.nz-2))
    return out
