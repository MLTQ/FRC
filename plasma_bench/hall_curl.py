"""Bounded cylindrical Faraday map and electric-field quadrature for Hall work."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.constants import mu_0 as MU0
from .toroidal_state import radius_squared


class HallCurl:
    """Geometry-only compatible work pairing; does not advance a fluid state."""
    def __init__(self, mesh, quadrature_order=3):
        if type(quadrature_order) is not int or quadrature_order < 2:
            raise ValueError('At least second-order Gauss electric quadrature required')
        self.mesh = mesh
        self.order = quadrature_order
        nr, nz = mesh.nr, mesh.nz
        x, w = leggauss(quadrature_order)
        shape = (nr, nz, quadrature_order, quadrature_order)
        self.r = np.broadcast_to(mesh.r[:, None, None, None]+.5*mesh.dr*x[None, None, :, None], shape)
        self.z = np.broadcast_to(mesh.z[None, :, None, None]+.5*mesh.dz*x[None, None, None, :], shape)
        self.weight = 2*np.pi*self.r*mesh.dr*mesh.dz/4*w[None, None, :, None]*w[None, None, None, :]
        i, j, _, _ = np.indices(shape)
        row = np.arange(np.prod(shape)).reshape(shape)
        az = (self.z-mesh.zf[j])/mesh.dz
        ar = (self.r-mesh.rf[i])/mesh.dr
        ax = (self.r**2-mesh.rf[i]**2)/(2*mesh.radial_measure[i])
        maps = []
        for component in range(3):
            rows, columns, data = [], [], []
            def add(col, value, mask=None):
                mask = np.ones(shape, bool) if mask is None else mask
                rows.append(row[mask]); columns.append(col[mask]); data.append(value[mask])
            if component == 0:
                add(i*nz+j, self.r*(1-az))
                add(i*nz+(j+1)%nz, self.r*az)
                size = nr*nz
            elif component == 1:
                for di, wr in ((0, 1-ar), (1, ar)):
                    ni = i+di
                    for dj, wz in ((0, 1-az), (1, az)):
                        add((ni-1)*nz+(j+dj)%nz, wr*wz, (ni > 0) & (ni < nr))
                size = (nr-1)*nz
            else:
                add(i*nz+j, 1-ax)
                add((i+1)*nz+j, ax, i+1 < nr)
                size = nr*nz
            maps.append(coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(columns))),
                                   shape=(row.size, size)).tocsr())
        self.maps = tuple(maps)
        self.mass = tuple((p.T @ p.multiply(self.weight.reshape(-1, 1))).tocsc() for p in self.maps)
        self.solves = tuple(splu(m).solve for m in self.mass)
        sizes = [p.shape[1] for p in self.maps]
        self.offsets = np.cumsum([0]+sizes)
        cells, nodes = nr*nz, (nr+1)*nz
        rows, cols, data = [], [], []
        for i in range(nr):
            for j in range(nz):
                cell = nodes+i*nz+j
                for col, value in ((i*nz+j, 1/mesh.dz), (i*nz+(j+1)%nz, -1/mesh.dz),
                                    (self.offsets[2]+i*nz+j, -1/mesh.radial_measure[i])):
                    rows.append(cell); cols.append(col); data.append(value)
                if i+1 < nr:
                    rows.append(cell); cols.append(self.offsets[2]+(i+1)*nz+j); data.append(1/mesh.radial_measure[i])
        for i in range(1, nr):
            for j in range(nz):
                rows.append(i*nz+j); cols.append(self.offsets[1]+(i-1)*nz+j); data.append(-mesh.rf[i])
        self.faraday = coo_matrix((data, (rows, cols)), shape=(nodes+cells, self.offsets[-1])).tocsr()

    def pack(self, psi, g):
        return np.concatenate((psi[:, :-1].ravel(), g.ravel()))

    def unpack(self, values):
        mesh = self.mesh
        nodes = (mesh.nr+1)*mesh.nz
        psi = np.asarray(values[:nodes]).reshape(mesh.nr+1, mesh.nz)
        return np.concatenate((psi, psi[:, :1]), axis=1), np.asarray(values[nodes:]).reshape(mesh.nr, mesh.nz)

    def split_electric(self, values):
        return tuple(values[self.offsets[j]:self.offsets[j+1]] for j in range(3))

    def sample_electric(self, values):
        parts = self.split_electric(values)
        return np.stack([p @ v for p, v in zip(self.maps, parts)], axis=-1).reshape(self.r.shape+(3,))

    def project_electric(self, values):
        values = np.asarray(values).reshape(-1, 3)
        return np.concatenate([solve(p.T @ (self.weight.ravel()*values[:, j]))
                               for j, (solve, p) in enumerate(zip(self.solves, self.maps))])

    def magnetic_gradient(self, psi, g):
        """Exact gradient of the existing cell-poloidal/annular-toroidal energy."""
        mesh = self.mesh
        br, bz = mesh.cell_fields(psi)
        wr = ((mesh.radial_centroid-mesh.rf[:-1])/mesh.dr)[:, None]
        ar = np.zeros((mesh.nr+1, mesh.nz))
        ar[:-1] += mesh.volumes*br*(1-wr)/MU0
        ar[1:] += mesh.volumes*br*wr/MU0
        az = np.zeros((mesh.nr, mesh.nz+1))
        az[:, :-1] += .5*mesh.volumes*bz/MU0
        az[:, 1:] += .5*mesh.volumes*bz/MU0
        radial = np.divide(ar, mesh.rf[:, None]*mesh.dz, out=np.zeros_like(ar), where=mesh.rf[:, None] > 0)
        axial = az/mesh.radial_measure[:, None]
        gradient = np.zeros_like(psi)
        gradient[:, :-1] += radial
        gradient[:, 1:] -= radial
        gradient[:-1] -= axial
        gradient[1:] += axial
        gradient[:, 0] += gradient[:, -1]
        return np.concatenate((gradient[:, :-1].ravel(), (mesh.volumes*radius_squared(mesh)*g/MU0).ravel()))

    def current(self, psi, g):
        dual = self.split_electric(-self.faraday.T @ self.magnetic_gradient(psi, g))
        coefficients = np.concatenate([solve(b) for solve, b in zip(self.solves, dual)])
        return self.sample_electric(coefficients)

    def magnetic_samples(self, psi, g):
        """Declared linear face-field interpolation, regular Btheta=r g."""
        mesh = self.mesh
        br, bz = mesh.face_fields(psi)
        ar = (self.r-mesh.rf[:-1, None, None, None])/mesh.dr
        az = (self.z-mesh.zf[None, :-1, None, None])/mesh.dz
        radial = (1-ar)*br[:-1, :, None, None]+ar*br[1:, :, None, None]
        axial = (1-az)*bz[:, :-1, None, None]+az*bz[:, 1:, None, None]
        return np.stack((radial, self.r*g[..., None, None], axial), axis=-1)
