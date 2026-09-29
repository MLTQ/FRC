"""Compatible Hall curl with prescribed zero toroidal field at the radial wall."""
import numpy as np
from scipy.sparse import coo_matrix, hstack
from scipy.sparse.linalg import splu
from .hall_compatible import CompatibleHallCurl


class TraceHallCurl(CompatibleHallCurl):
    """Fixed psi, Btheta(R)=0; Ez(R) is free, not a conducting electric trace."""
    wall_btheta = 0.

    def __init__(self, mesh, quadrature_order=5):
        super().__init__(mesh, quadrature_order)
        nr, nz = mesh.nr, mesh.nz
        i, j, _, _ = np.indices(self.r.shape)
        row = np.arange(self.r.size).reshape(self.r.shape)
        x = (self.r**2-mesh.rf[i]**2)/(2*mesh.radial_measure[i])
        axial = coo_matrix((np.concatenate(((1-x).ravel(), x.ravel())),
                            (np.tile(row.ravel(), 2), np.concatenate(((i*nz+j).ravel(), ((i+1)*nz+j).ravel())))),
                           shape=(self.r.size, (nr+1)*nz)).tocsr()
        mass = (axial.T @ axial.multiply(self.weight.reshape(-1, 1))).tocsc()
        self.maps = (*self.maps[:2], axial)
        self.mass = (*self.mass[:2], mass)
        self.solves = (*self.solves[:2], splu(mass).solve)
        nodes = (nr+1)*nz
        row = nodes+(nr-1)*nz+np.arange(nz)
        boundary = coo_matrix((np.full(nz, 1/mesh.radial_measure[-1]), (row, np.arange(nz))),
                              shape=(nodes+nr*nz, nz)).tocsr()
        self.faraday = hstack((self.faraday, boundary), format='csr')
        self.offsets = np.cumsum([0]+[p.shape[1] for p in self.maps])

    def scalar_samples(self, values):
        """Positive Q1 sampling of a full periodic dual scalar at Hall points."""
        values = np.asarray(values, float)
        mesh = self.mesh
        if values.shape != (mesh.nr+1, mesh.nz) or not np.isfinite(values).all():
            raise ValueError('Finite full dual scalar required')
        x = (self.r**2-mesh.rf[:-1, None, None, None]**2)/(2*mesh.radial_measure[:, None, None, None])
        lo = (1-self.az)*values[:-1, :, None, None]+self.az*np.roll(values[:-1], -1, axis=1)[..., None, None]
        hi = (1-self.az)*values[1:, :, None, None]+self.az*np.roll(values[1:], -1, axis=1)[..., None, None]
        return (1-x)*lo+x*hi
