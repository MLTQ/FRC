"""Galerkin magnetic space with the intended outer toroidal trace imposed."""
from copy import copy
import numpy as np
from scipy.sparse import block_diag, eye, kron, csr_matrix
from scipy.sparse.linalg import splu
from .compatible_geometry import CompatibleGeometry, gram
from .compatible_polynomials import bases


class StrongWallGeometry(CompatibleGeometry):
    def __init__(self, original):
        self.__dict__.update(copy(original.__dict__))
        self.original = original
        mesh = self.mesh
        self.wall_trace = bases(.5*mesh.rf**2, np.array([.5*mesh.radius**2]), self.degree)[2].toarray()[0]
        radial_mass = gram(self.radial_e, self.wr)
        dual = splu(radial_mass).solve(self.wall_trace)
        self.radial_projection = np.eye(mesh.nr)-np.outer(dual, self.wall_trace)/(self.wall_trace@dual)
        nodes = (mesh.nr+1)*mesh.nz
        self.projection = block_diag((eye(nodes), kron(csr_matrix(self.radial_projection), eye(mesh.nz))), format='csr')
        self.faraday = (self.projection@original.faraday).tocsr()

    def project_curl(self, value):
        value = np.asarray(value, float); count = (self.mesh.nr+1)*self.mesh.nz
        return np.r_[value[:count], (self.radial_projection@value[count:].reshape(self.mesh.nr, self.mesh.nz)).ravel()]

    def wall_residual(self, g):
        g = np.asarray(g, float)
        if g.shape != (self.mesh.nr, self.mesh.nz): raise ValueError('Complete toroidal coefficients required')
        scale = np.linalg.norm(self.wall_trace)*np.linalg.norm(g)
        return float(np.linalg.norm(self.wall_trace@g)/max(scale, np.finfo(float).tiny))

    def pack(self, psi, g):
        value = super().pack(psi, g)
        if self.wall_residual(g) > 256*np.finfo(float).eps:
            raise ValueError('Boundary-conforming toroidal state required; no state repair')
        return value
