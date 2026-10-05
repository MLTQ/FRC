"""Fixed physical scaling and active coordinates for coupled implicit solves."""
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from .regular_plasma_step import field_from_packed


class PlasmaCoordinates:
    def __init__(self, system, state, field):
        self.system = system; self.fluid = d = system.fluid; self.geometry = system.geometry
        f = d.fields(state); field.packed_increment(self.geometry)
        self.bias = field.axial_bias
        self.density = float(np.sum(d.volume*state[..., 0])/np.sum(d.volume))
        self.pressures = np.sum(d.volume[..., None]*f['pressure'], axis=(0, 1, 2))/np.sum(d.volume)
        self.speed = float(np.sqrt(np.sum(d.thermal.gammas*self.pressures)/(d.mass*self.density)))
        self.axial_speed = float(np.sum(d.volume*state[..., 3])/np.sum(d.volume*state[..., 0]))
        active = np.ones(d.shape, bool); active[0, 0, :, 1:3] = False
        self.fluid_indices = np.flatnonzero(active)
        mesh = self.geometry.mesh; count = (mesh.nr+1)*mesh.nz
        self.magnetic_indices = np.r_[np.arange(mesh.nz, mesh.nr*mesh.nz),
                                      np.arange(count, count+mesh.nr*mesh.nz)]
        scales = self.density*np.array([1., self.speed/d.radius, self.speed/d.radius, self.speed, 1., 1.])
        self.fluid_scale = np.broadcast_to(scales, d.shape).ravel()[self.fluid_indices]
        bscale = np.sqrt(MU0*d.mass*self.density)*self.speed
        self.magnetic_scale = np.r_[np.full(count, bscale*d.radius**2),
                                     np.full(mesh.nr*mesh.nz, bscale/d.radius)][self.magnetic_indices]
        self.scale = np.r_[self.fluid_scale, self.magnetic_scale]
        self.nfluid = len(self.fluid_indices); self.size = len(self.scale)
        if not np.isfinite(self.scale).all() or np.any(self.scale <= 0):
            raise ValueError('Positive finite physical coordinate scales required')

    def pack(self, fluid, magnetic):
        if np.shape(fluid) != self.fluid.shape or np.shape(magnetic) != (self.geometry.hessian.shape[0],):
            raise ValueError('Complete fluid and packed magnetic shapes required')
        return np.r_[np.asarray(fluid).ravel()[self.fluid_indices],
                     np.asarray(magnetic)[self.magnetic_indices]]/self.scale

    def unpack(self, value):
        value = np.asarray(value, float)
        if value.shape != (self.size,) or not np.isfinite(value).all():
            raise ValueError('Finite complete active coordinate vector required')
        fluid = np.zeros(self.fluid.shape); magnetic = np.zeros(self.geometry.hessian.shape[0])
        fluid.ravel()[self.fluid_indices] = value[:self.nfluid]*self.fluid_scale
        magnetic[self.magnetic_indices] = value[self.nfluid:]*self.magnetic_scale
        return fluid, magnetic

    def background(self):
        d = self.fluid; u = np.zeros(d.shape); u[..., 0] = self.density
        u[..., 3] = self.density*self.axial_speed
        entropy = np.log(self.pressures/(E*self.density))-(d.thermal.gammas-1)*np.log(self.density/d.thermal.n0)
        u[..., 4:] = self.density*entropy
        return u

    def field(self, packed):
        return field_from_packed(self.geometry, packed, self.bias)
