"""Finite axial polynomial maps and physical derivatives on the actual fluid."""
import numpy as np
from scipy.sparse import csr_matrix
from .regular_field_map import RegularFieldMap
from .physical_derivatives import fluid_matrices,RegularDifferential
from .regular_species import RegularSpecies
from .variable_compatible_current import weak_vector


def axial_matrix(fluid,points,*,derivative=False):
    nodes=2*fluid.z.ravel()/fluid.length;points=2*np.asarray(points)/fluid.length
    v=np.polynomial.legendre.legvander(nodes,fluid.nz-1)
    values=np.polynomial.legendre.legvander(points,fluid.nz-1)
    if derivative:
        values=np.column_stack([np.polynomial.legendre.Legendre.basis(k).deriv()(points) for k in range(fluid.nz)])*2/fluid.length
    return csr_matrix(np.linalg.solve(v.T,values.T).T)


class FiniteCylinderMap(RegularFieldMap):
    def __init__(self,fluid,geometry,*,side=None):
        m=geometry.mesh
        if not np.isclose(m.radius,fluid.radius,rtol=1e-14,atol=0) or not np.isclose(2*m.half_length,fluid.length,rtol=1e-14,atol=0) or m.nr%fluid.elements:
            raise ValueError('Matching finite cylinder and subdivided radial fluid elements required')
        self.fluid=fluid;self.geometry=geometry;self.shape=fluid.shape[:-1]
        self.radial=fluid_matrices(fluid,geometry.radial_points,side)
        self.radial_velocity=fluid_matrices(fluid,geometry.radial_points,side,transverse=True)
        self.axial=axial_matrix(fluid,geometry.axial_points)


class FiniteFluidDifferential(RegularDifferential):
    def __init__(self,mapping):
        self.mapping=mapping;self.geometry=mapping.geometry;d=mapping.fluid
        self.dr=fluid_matrices(d,self.geometry.radial_points,derivative=True)
        self.drv=fluid_matrices(d,self.geometry.radial_points,derivative=True,transverse=True)
        self.dz=axial_matrix(d,self.geometry.axial_points,derivative=True)


class FiniteSpecies(RegularSpecies):
    def __init__(self,model):
        super().__init__(model)
        self.potential.shape=(model.mesh.nr+1,model.mesh.nz+1)

    def gradient_adjoint(self,vector):
        g=self.geometry;m=self.model.mesh;er,_,ez=g.split_electric(weak_vector(g,vector))
        value=g.Dx.T@er.reshape(m.nr,m.nz+1)
        value+=ez.reshape(m.nr+1,m.nz)@g.Dz
        return self.potential.adjoint(np.asarray(value))
