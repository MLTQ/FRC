"""Separable physical wall coordinates without a dense full-grid null space."""
import numpy as np
from scipy.linalg import null_space
from .plasma_coordinates import PlasmaCoordinates


class ModalWallCoordinates:
    def __init__(self, system, state, field):
        self.full = full = PlasmaCoordinates(system.original_system, state, field)
        self.system = system; self.fluid = d = system.fluid; self.geometry = geo = system.geometry
        if d.nz != geo.mesh.nz: raise ValueError('Matching axial fluid/field counts required for modal inverse')
        d.fields(state); field.packed_increment(geo)
        for key in ('bias', 'density', 'pressures', 'speed', 'axial_speed'):
            setattr(self, key, getattr(full, key))
        self.nz = d.nz; points = d.elements*(d.degree+1)
        active = np.ones((points, 6), bool); active[0, 1:3] = False; active[-1, 1] = False
        self.fluid_rows = np.flatnonzero(active)
        scales = self.density*np.array([1., self.speed/d.radius, self.speed/d.radius, self.speed, 1., 1.])
        self.fluid_row_scale = np.tile(scales, points)[self.fluid_rows]
        self.radial_basis = null_space(geo.wall_trace[None, :])
        nr = geo.mesh.nr; self.magnetic_scalars = np.array([full.magnetic_scale[0], full.magnetic_scale[-1]])
        self.fluid_row_count = len(self.fluid_rows); self.magnetic_row_count = 2*(nr-1)
        self.nfluid = self.fluid_row_count*self.nz
        self.size = (self.fluid_row_count+self.magnetic_row_count)*self.nz

    def pack(self, fluid, magnetic):
        if np.shape(fluid) != self.fluid.shape or np.shape(magnetic) != (self.geometry.hessian.shape[0],):
            raise ValueError('Complete fluid/magnetic arrays required')
        f = np.asarray(fluid).reshape(-1, self.nz, 6).transpose(0, 2, 1).reshape(-1, self.nz)
        psi, g = self.geometry.unpack(np.asarray(magnetic))
        return np.r_[(f[self.fluid_rows]/self.fluid_row_scale[:, None]).ravel(),
            (psi[1:-1, :-1]/self.magnetic_scalars[0]).ravel(),
            (self.radial_basis.T@g/self.magnetic_scalars[1]).ravel()]

    def unpack(self, value):
        value = np.asarray(value, float)
        if value.shape != (self.size,) or not np.isfinite(value).all():
            raise ValueError('Finite full modal-wall coordinate vector required')
        rows = value.reshape(-1, self.nz); nf = self.fluid_row_count; nr = self.geometry.mesh.nr
        f = np.zeros((self.fluid.elements*(self.fluid.degree+1)*6, self.nz))
        f[self.fluid_rows] = rows[:nf]*self.fluid_row_scale[:, None]
        fluid = f.reshape(-1, 6, self.nz).transpose(0, 2, 1).reshape(self.fluid.shape)
        psi = np.zeros((nr+1, self.nz+1)); psi[1:-1, :-1] = rows[nf:nf+nr-1]*self.magnetic_scalars[0]
        psi[:, -1] = psi[:, 0]
        g = self.radial_basis@rows[nf+nr-1:]*self.magnetic_scalars[1]
        return fluid, self.geometry.pack(psi, g)

    def background(self): return self.full.background()

    def field(self, packed):
        result = self.full.field(packed); result.packed_increment(self.geometry)
        return result

    def fourier(self, value):
        return np.fft.fft(np.asarray(value).reshape(-1, self.nz), axis=1)

    def physical(self, modes):
        return np.fft.ifft(modes, axis=1).ravel()
