"""Physical coordinates parameterizing every boundary-conforming stage."""
import numpy as np
from scipy.linalg import null_space
from .plasma_coordinates import PlasmaCoordinates


class StrongWallCoordinates:
    def __init__(self, system, state, field):
        full = PlasmaCoordinates(system.original_system, state, field)
        self.__dict__.update(full.__dict__)
        self.full = full; self.system = system; self.fluid = system.fluid
        self.geometry = geometry = system.geometry
        system.fluid.fields(state); field.packed_increment(geometry)
        wall = np.zeros(self.fluid.shape, bool); wall[-1, -1, :, 1] = True
        self.fluid_keep = np.flatnonzero(~wall.ravel()[full.fluid_indices])
        count = (geometry.mesh.nr+1)*geometry.mesh.nz
        constraint = np.zeros((geometry.mesh.nz, geometry.hessian.shape[0]))
        constraint[:, count:] = np.kron(geometry.wall_trace[None, :], np.eye(geometry.mesh.nz))
        self.constraint = constraint[:, full.magnetic_indices]*full.magnetic_scale
        self.magnetic_basis = null_space(self.constraint)
        self.nfluid = len(self.fluid_keep)
        self.size = self.nfluid+self.magnetic_basis.shape[1]

    def pack(self, fluid, magnetic):
        value = self.full.pack(fluid, magnetic)
        return np.r_[value[:self.full.nfluid][self.fluid_keep],
                     self.magnetic_basis.T@value[self.full.nfluid:]]

    def unpack(self, value):
        value = np.asarray(value, float)
        if value.shape != (self.size,) or not np.isfinite(value).all():
            raise ValueError('Finite complete wall-conforming coordinate vector required')
        full = np.zeros(self.full.size)
        full[self.fluid_keep] = value[:self.nfluid]
        full[self.full.nfluid:] = self.magnetic_basis@value[self.nfluid:]
        return self.full.unpack(full)

    def field(self, packed):
        result = self.full.field(packed)
        result.packed_increment(self.geometry)
        return result

    def background(self):
        return self.full.background()
