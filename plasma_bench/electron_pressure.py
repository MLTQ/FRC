"""Paired adiabatic electron entropy drift and pressure-driven toroidal induction."""
import numpy as np
from scipy.constants import elementary_charge as E
from .dual_current import DualCurrent
from .thermal_means import entropy_mean


class ElectronPressure:
    """Fixed positive density split component, not a complete two-fluid solver."""
    def __init__(self, mesh, number_density, *, gamma=5/3, density_scale=1e20, wall_btheta=0.):
        self.geometry = DualCurrent(mesh)
        self.density = self.geometry.scalar(number_density, 'electron density').copy()
        if np.any(self.density <= 0) or not np.isfinite([gamma, density_scale, wall_btheta]).all() or gamma <= 1 or density_scale <= 0:
            raise ValueError('Positive density/scale, gamma > 1 and finite constant wall Btheta required')
        if np.ptp(self.density[-1]) > 1e-13*self.density[-1].max():
            raise ValueError('This closed pressure component requires constant wall density')
        self.gamma, self.density_scale, self.wall_btheta = gamma, density_scale, float(wall_btheta)
        self.a = (self.density/density_scale)**(gamma-1)
        self.energy_weight = self.a/(gamma-1)
        for values in (self.density, self.a, self.energy_weight):
            values.setflags(write=False)

    def from_pressure(self, pressure):
        pressure = self.geometry.scalar(pressure, 'electron pressure')
        if np.any(pressure <= 0):
            raise ValueError('Positive electron pressure required; no floor applied')
        return pressure/self.a

    def pressure(self, entropy_density):
        return self.a*self.geometry.scalar(entropy_density, 'entropy proxy density')

    def rhs(self, g, entropy_density):
        """U=n w, p=a U; w is a scaled advected polytropic entropy proxy."""
        geo = self.geometry
        u = geo.scalar(entropy_density, 'entropy proxy density')
        if np.any(u <= 0):
            raise ValueError('Nonpositive electron thermal state; no repair applied')
        w = u/self.density
        wr, wz = entropy_mean(w[1:], w[:-1]), entropy_mean(w, np.roll(w, -1, axis=1))
        tr, tz = geo.gradient(self.a*w)
        ar, az = geo.gradient(self.energy_weight)
        # -grad(p/n)/e - w grad(a/(gamma-1))/e equals -grad(p)/(n e).
        er, ez = -(tr+wr*ar)/E, -(tz+wz*az)/E
        # C grad(T)=0 analytically. Subtract one local constant w from the
        # closed A circulation too, avoiding cancellation of large gradients.
        da_r = np.diff(self.energy_weight, axis=0)
        da_z = np.roll(self.energy_weight, -1, axis=1)-self.energy_weight
        anchor = wz[:-1]
        circulation = ((wz[1:]-anchor)*da_z[1:]
                       -(np.roll(wr, -1, axis=1)-anchor)*np.roll(da_r, -1, axis=1)
                       +(wr-anchor)*da_r)
        dg = -circulation/(E*geo.dx[:, None]*geo.mesh.dz)
        ir, iz = geo.current_flux(g, self.wall_btheta)
        radial_flux = np.zeros_like(ir)
        radial_flux[1:-1] = ir[1:-1]*wr/E
        du = geo.divergence(radial_flux, iz*wz/E)
        return dg, du, er, ez

    def diagnostics(self, g, entropy_density):
        geo = self.geometry
        u = geo.scalar(entropy_density)
        p = self.pressure(u)
        if np.any(p <= 0):
            raise ValueError('Positive electron pressure required for diagnostics')
        magnetic = float(.5*np.sum(geo.magnetic_mass*g*g))
        thermal = float(np.sum(geo.volumes*self.energy_weight*u))
        w = u/self.density
        # Fixed SI reference affects the constant entropy gauge, not its change.
        entropy = float(np.sum(geo.volumes*self.density*np.log(w/E)))
        ir, iz = geo.current_flux(g, self.wall_btheta)
        return dict(magnetic_energy_j=magnetic, electron_internal_energy_j=thermal,
                    total_energy_j=magnetic+thermal, entropy_proxy_integral_j=float(np.sum(geo.volumes*u)),
                    entropy_log_particle_units=entropy, min_electron_pressure_pa=float(p.min()),
                    max_charge_divergence_a_m3=float(abs(geo.divergence(ir, iz)).max()))
