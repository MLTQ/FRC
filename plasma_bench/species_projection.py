"""Physical scalar L2 transfer and transpose between regular and field spaces."""
import numpy as np
from scipy.sparse.linalg import splu
from .compatible_geometry import gram


class SpeciesProjection:
    """T=(N.T W N)^-1 N.T W S; no nodal temperature or density replacement."""
    def __init__(self, model):
        self.mapping = model.mapping
        geo = model.geometry
        self.radial = splu(gram(geo.Nx, geo.wr)).solve(
            (geo.Nx.T @ self.mapping.radial.multiply(geo.wr[:, None])).toarray())
        self.axial = splu(gram(geo.Nz, geo.wz)).solve(
            (geo.Nz.T @ self.mapping.axial.multiply(geo.wz[:, None])).toarray())
        self.shape = (model.mesh.nr+1, model.mesh.nz)

    def apply(self, value):
        value = self.mapping.validate(value)
        return self.radial @ value.reshape(-1, self.mapping.fluid.nz) @ self.axial.T

    def adjoint(self, value):
        if np.shape(value) != self.shape or not np.isfinite(value).all():
            raise ValueError('Finite compatible scalar covector required')
        return (self.radial.T @ value @ self.axial).reshape(self.mapping.shape)
