"""Original tensor-product geometry with separate one-dimensional mass solves."""
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.sparse import block_diag, diags, kron
from scipy.sparse.linalg import splu
from .compatible_geometry import CompatibleGeometry, gram
from .compatible_polynomials import bases


class TensorMassInverse:
    def __init__(self, radial, axial):
        self.shape = (radial.shape[0], axial.shape[0])
        self.radial = splu(radial); self.axial = splu(axial)

    def solve(self, value):
        value = np.asarray(value)
        if value.shape != (self.shape[0]*self.shape[1],) or np.iscomplexobj(value) or not np.isfinite(value).all():
            raise ValueError('Complete finite real mass-dual vector required')
        first = self.radial.solve(value.reshape(self.shape))
        return self.axial.solve(first.T).T.ravel()


class TensorMassGeometry(CompatibleGeometry):
    def __init__(self, mesh, *, degree=5, quadrature_order=10):
        if type(quadrature_order) is not int or quadrature_order < degree+1:
            raise ValueError('Physical quadrature must resolve the selected polynomial degree')
        self.mesh = mesh; self.degree = degree; self.order = quadrature_order
        q, w = np.polynomial.legendre.leggauss(quadrature_order)
        self.radial_points = (mesh.r[:, None]+.5*mesh.dr*q).ravel()
        self.axial_points = (mesh.z[:, None]+.5*mesh.dz*q).ravel()
        self.wr = (np.pi*mesh.dr*self.radial_points.reshape(mesh.nr, -1)*w).ravel()
        self.wz = np.broadcast_to(.5*mesh.dz*w, (mesh.nz, len(q))).ravel()
        self.Nx, self.dNx, self.Ex, self.Dx = bases(.5*mesh.rf**2, .5*self.radial_points**2, degree)
        self.Nz, self.dNz, self.Ez, self.Dz = bases(mesh.zf, self.axial_points, degree, periodic=True)
        radius = self.radial_points
        self.radial_e = self.Ex.multiply(radius[:, None]).tocsr()
        self.theta_e = (self.Nx[:, 1:-1].multiply(1/radius[:, None])@diags(mesh.rf[1:-1])).tocsr()
        self.radial_b = self.Nx.multiply(1/radius[:, None]).tocsr()
        self.radial_b = self.radial_b@diags(np.r_[0., np.ones(mesh.nr)])
        self.pairs = ((self.radial_e, self.Nz), (self.theta_e, self.Nz), (self.Nx, self.Ez))
        factors = tuple((gram(a, self.wr), gram(b, self.wz)) for a, b in self.pairs)
        self.mass = tuple(kron(a, b, format='csc') for a, b in factors)
        self.tensor_mass_factors = tuple(TensorMassInverse(a, b) for a, b in factors)
        self.solves = tuple(factor.solve for factor in self.tensor_mass_factors)
        self.offsets = np.cumsum([0]+[m.shape[0] for m in self.mass])
        hpsi = kron(gram(self.radial_b, self.wr), gram(self.dNz, self.wz))+kron(gram(self.dNx, self.wr), gram(self.Nz, self.wz))
        hg = kron(gram(self.radial_e, self.wr), gram(self.Ez, self.wz))
        self.hessian = block_diag((hpsi, hg), format='csc')/MU0
        self.faraday = self._faraday()
        self.r = np.broadcast_to(radius.reshape(mesh.nr, 1, len(q), 1), (mesh.nr, mesh.nz, len(q), len(q)))
        self.z = np.broadcast_to(self.axial_points.reshape(1, mesh.nz, 1, len(q)), self.r.shape)
        self.weight = self._layout(self.wr[:, None]*self.wz[None])
