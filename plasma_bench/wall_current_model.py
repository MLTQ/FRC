"""Physical wall current and exact elimination of its coupled inertia block."""
from copy import copy
from dataclasses import dataclass
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse.linalg import LinearOperator
from .regular_current import RegularCurrent
from .regular_field_map import RegularFieldMap


@dataclass(frozen=True)
class WallCurrentField:
    magnetic: object
    wall_current: np.ndarray

    @property
    def axial_bias(self): return self.magnetic.axial_bias
    def packed_increment(self, geometry): return self.magnetic.packed_increment(geometry)
    def samples(self, geometry): return self.magnetic.samples(geometry)


class WallCurrentModel(RegularCurrent):
    def __init__(self, mapping, number_density, ion_mass, wall_current):
        super().__init__(mapping,number_density,ion_mass)
        self.wall_current=np.asarray(wall_current,float)
        if self.wall_current.shape!=(self.mesh.nz,) or not np.isfinite(self.wall_current).all():
            raise ValueError('Finite complete evolved wall current required')
        base_mapping=RegularFieldMap(mapping.fluid,self.geometry.base_geometry)
        self.base=RegularCurrent(base_mapping,number_density,ion_mass)

    def current(self, magnetic):
        return self.geometry.embed(self.base.current(magnetic),self.wall_current)

    def inertia_blocks(self):
        """Actual I0chi and Ichichi, including changing/nonuniform density."""
        g=self.geometry;b=g.base_geometry;chi=g.theta_e[:,-1].toarray().ravel()
        coefficient=g._unlayout(self.inertia_coefficient)
        radial=b.theta_e.T@((g.wr*chi)[:,None]*coefficient)
        z=g.Nz.toarray();wz=g.wz
        cross=np.zeros((b.offsets[-1],self.mesh.nz))
        for k,row in enumerate(radial):
            lo=b.offsets[1]+k*self.mesh.nz
            cross[lo:lo+self.mesh.nz]=z.T@((wz*row)[:,None]*z)
        wall=z.T@((wz*((g.wr*chi*chi)@coefficient))[:,None]*z)
        return cross,cho_factor(wall,lower=True)

    def reduced_induction(self, forcing):
        """Schur elimination of the physical wall-current time derivative."""
        g=self.geometry;cross,factor=self.inertia_blocks()
        wall_force=forcing[g.wall_indices]
        solve=lambda x:cho_solve(factor,x)
        reduced=copy(self.base);original=self.base.inertia
        def action(value):return original@value-cross@solve(cross.T@value)
        reduced.inertia=LinearOperator(original.shape,matvec=action,rmatvec=action,dtype=float)
        right=forcing[g.old_indices]-cross@solve(wall_force)
        def reconstruct(old_rate):
            wall_rate=-solve(wall_force+cross.T@old_rate)
            return g.embed(old_rate,wall_rate),wall_rate
        return reduced,right,reconstruct
