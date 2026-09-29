"""Cylindrical Hall electric and magnetic spaces with a commuting Faraday curl."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.constants import mu_0 as MU0
from .hall_curl import HallCurl
from .compatible_magnetic import energy_density, energy_rate, gradient
from .toroidal_state import radius_squared


class CompatibleHallCurl(HallCurl):
    """Use psi/r electric shapes and the explicitly compatible magnetic energy."""
    def __init__(self, mesh, quadrature_order=5):
        super().__init__(mesh, quadrature_order)
        i, j, _, _ = np.indices(self.r.shape)
        row = np.arange(self.r.size).reshape(self.r.shape)
        ax = (self.r**2-mesh.rf[i]**2)/(2*mesh.radial_measure[i])
        self.az = (self.z-mesh.zf[j])/mesh.dz
        self.left_radial = mesh.rf[i]/self.r*(1-ax)
        self.right_radial = mesh.rf[i+1]/self.r*ax
        rows, columns, values = [], [], []
        for di, wr in ((0, self.left_radial), (1, self.right_radial)):
            ni = i+di
            valid = (ni > 0) & (ni < mesh.nr)
            for dj, wz in ((0, 1-self.az), (1, self.az)):
                rows.append(row[valid])
                columns.append(((ni-1)*mesh.nz+(j+dj)%mesh.nz)[valid])
                values.append((wr*wz)[valid])
        theta = coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))),
                           shape=self.maps[1].shape).tocsr()
        mass = (theta.T @ theta.multiply(self.weight.reshape(-1, 1))).tocsc()
        self.maps = (self.maps[0], theta, self.maps[2])
        self.mass = (self.mass[0], mass, self.mass[2])
        self.solves = (self.solves[0], splu(mass).solve, self.solves[2])

    def magnetic_samples(self, psi, g):
        br, bz = self.mesh.face_fields(psi)
        radial = self.left_radial*br[:-1, :, None, None]+self.right_radial*br[1:, :, None, None]
        axial = (1-self.az)*bz[:, :-1, None, None]+self.az*bz[:, 1:, None, None]
        return np.stack((radial, self.r*g[..., None, None], axial), axis=-1)

    def magnetic_gradient(self, psi, g):
        mesh = self.mesh
        return np.concatenate((gradient(mesh, psi)[:, :-1].ravel(),
                               (mesh.volumes*radius_squared(mesh)*g/MU0).ravel()))

    def energy_density(self, psi, g):
        return energy_density(self.mesh, psi)+radius_squared(self.mesh)*g*g/(2*MU0)

    def energy_rate(self, psi, g, psi_rate, g_rate):
        return energy_rate(self.mesh, psi, psi_rate)+radius_squared(self.mesh)*g*g_rate/MU0
