"""Independent cylindrical pressure-battery and electron drift references."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0


@dataclass(frozen=True)
class ElectronReference:
    radius: float = .035
    length: float = .24
    density: float = 1e20
    temperature_ev: float = 5.
    density_contrast: float = .4
    temperature_contrast: float = .3
    toroidal_slope: float = .05
    gamma: float = 5/3

    def fields(self, r, z):
        r, z = np.broadcast_arrays(r, z)
        s, k = (r/self.radius)**2, 2*np.pi/self.length
        f, fs = (1-s)**2, -2*(1-s)
        h, hs = (1-s)**3, -3*(1-s)**2
        cz, sz = np.cos(k*z), np.sin(k*z)
        ln_r_over_r = self.density_contrast*2/self.radius**2*fs*cz
        ln_z = -self.density_contrast*k*f*sz
        n = self.density*np.exp(self.density_contrast*f*cz)
        t = self.temperature_ev*E*(1+self.temperature_contrast*h*sz)
        t_r_over_r = self.temperature_ev*E*self.temperature_contrast*2/self.radius**2*hs*sz
        t_z = self.temperature_ev*E*self.temperature_contrast*k*h*cz
        pressure = n*t
        er_over_r = -(t_r_over_r+t*ln_r_over_r)/E
        ez = -(t_z+t*ln_z)/E
        # Direct analytic curl of -grad(p)/(n e), including its regular-axis limit.
        g_rate = (ln_r_over_r*t_z-ln_z*t_r_over_r)/E
        a = (n/self.density)**(self.gamma-1)
        w = t/a
        w_r = r*(t_r_over_r-(self.gamma-1)*t*ln_r_over_r)/a
        w_z = (t_z-(self.gamma-1)*t*ln_z)/a
        phase = k*z+.37
        g = self.toroidal_slope*h*np.cos(phase)
        jr = self.toroidal_slope*k*r*h*np.sin(phase)/MU0
        jz = 2*self.toroidal_slope*(h+s*hs)*np.cos(phase)/MU0
        return dict(density=n, pressure=pressure, temperature_j=t, entropy_density=n*w,
                    g=g, g_rate=g_rate, entropy_rate=(jr*w_r+jz*w_z)/E,
                    jr=jr, jz=jz, er_over_r=er_over_r, ez=ez,
                    stream=2*np.pi*r*r*g/MU0)

    def cell_average(self, mesh, key, order=10):
        nodes, weights = np.polynomial.legendre.leggauss(order)
        lo, hi = .5*mesh.rf[:-1]**2, .5*mesh.rf[1:]**2
        result = np.zeros((mesh.nr, mesh.nz))
        for x, wx in zip(nodes, weights):
            r = np.sqrt(lo+hi+x*(hi-lo))[:, None]
            for z, wz in zip(nodes, weights):
                result += .25*wx*wz*self.fields(r, mesh.z[None]+.5*mesh.dz*z)[key]
        return result


def drift_pressure(geometry, time, *, density=1e20, slope=.01, pressure=80., amplitude=.3):
    velocity = -2*slope/(MU0*density*E)
    k = 2*np.pi/(2*geometry.mesh.half_length)
    return pressure*(1+amplitude*np.cos(k*(geometry.z-velocity*time)))
