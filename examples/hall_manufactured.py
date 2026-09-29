"""Exact nonlinear Hall induction and local Poynting identity on a filled cylinder."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0


@dataclass(frozen=True)
class HallManufactured:
    radius: float = .035
    length: float = .24
    number_density: float = 1e20
    density_contrast: float = .4
    bias: float = .005
    poloidal_amplitude: float = .002  # T
    toroidal_slope: float = .05      # T/m
    toroidal_mode: int = 2

    def __post_init__(self):
        if not np.isfinite([self.radius, self.length, self.number_density, self.density_contrast,
                            self.bias, self.poloidal_amplitude, self.toroidal_slope]).all():
            raise ValueError('Finite manufactured-field inputs required')
        if min(self.radius, self.length, self.number_density) <= 0 or abs(self.density_contrast) >= 1:
            raise ValueError('Positive geometry/density and |density contrast| < 1 required')
        if type(self.toroidal_mode) is not int or self.toroidal_mode < 1:
            raise ValueError('Positive integer toroidal axial mode required')

    def fields(self, r, z):
        """Analytic field/current and first derivatives, including the regular axis."""
        r, z = np.broadcast_arrays(r, z)
        x, k = r*r/self.radius**2, 2*np.pi/self.length
        if np.any(r < 0) or np.any(r > self.radius) or not np.isfinite(r).all() or not np.isfinite(z).all():
            raise ValueError('Finite coordinates inside the cylinder required')
        f, fp, fpp, fppp = (1-x)**3, -3*(1-x)**2, 6*(1-x), -6.
        a, c, ell = self.poloidal_amplitude, self.toroidal_slope, self.toroidal_mode*k
        sz, cz = np.sin(k*z), np.cos(k*z)
        st, ct = np.sin(ell*z+.37), np.cos(ell*z+.37)
        h, h2 = f+x*fp, 2*fp+x*fpp
        jt = a*k*k*r*f-4*a*r/self.radius**2*h2
        jt_r = a*k*k*(f+2*x*fp)-4*a/self.radius**2*(2*fp+7*x*fpp+2*x*x*fppp)
        magnetic = np.stack((a*k*r*f*sz, c*r*f*ct, self.bias+2*a*h*cz), axis=-1)
        magnetic_r = np.stack((a*k*(f+2*x*fp)*sz, c*(f+2*x*fp)*ct, 4*a*r/self.radius**2*h2*cz), axis=-1)
        magnetic_z = np.stack((a*k*k*r*f*cz, -c*ell*r*f*st, -2*a*k*h*sz), axis=-1)
        current = np.stack((c*ell*r*f*st, jt*cz, 2*c*h*ct), axis=-1)/MU0
        current_r = np.stack((c*ell*(f+2*x*fp)*st, jt_r*cz, 4*c*r/self.radius**2*h2*ct), axis=-1)/MU0
        current_z = np.stack((c*ell*ell*r*f*ct, -k*jt*sz, -2*c*ell*h*st), axis=-1)/MU0
        cn, sn = np.cos(k*z+.61), np.sin(k*z+.61)
        n = self.number_density*(1+self.density_contrast*f*cn)
        n_r = self.number_density*self.density_contrast*2*r/self.radius**2*fp*cn
        n_z = -self.number_density*self.density_contrast*k*f*sn
        electric = np.cross(current, magnetic)/(n[..., None]*E)
        electric_r = (np.cross(current_r, magnetic)+np.cross(current, magnetic_r))/(n[..., None]*E)-electric*(n_r/n)[..., None]
        electric_z = (np.cross(current_z, magnetic)+np.cross(current, magnetic_z))/(n[..., None]*E)-electric*(n_z/n)[..., None]
        et_over_r = np.divide(electric[..., 1], r, out=np.asarray(electric_r[..., 1]).copy(), where=r > 0)
        induction = np.stack((electric_z[..., 1], electric_r[..., 2]-electric_z[..., 0],
                              -electric_r[..., 1]-et_over_r), axis=-1)
        poynting = np.cross(electric, magnetic)/MU0
        poynting_r = (np.cross(electric_r, magnetic)+np.cross(electric, magnetic_r))/MU0
        poynting_z = (np.cross(electric_z, magnetic)+np.cross(electric, magnetic_z))/MU0
        sr_over_r = np.divide(poynting[..., 0], r, out=np.asarray(poynting_r[..., 0]).copy(), where=r > 0)
        return dict(magnetic=magnetic, magnetic_r=magnetic_r, magnetic_z=magnetic_z,
                    current=current, current_r=current_r, current_z=current_z,
                    number_density=n, electric=electric, electric_r=electric_r, electric_z=electric_z,
                    induction=induction, psi=.5*self.bias*r*r+a*r*r*f*cz,
                    psi_rate=-r*electric[..., 1], poynting=poynting,
                    poynting_divergence=sr_over_r+poynting_r[..., 0]+poynting_z[..., 2])
