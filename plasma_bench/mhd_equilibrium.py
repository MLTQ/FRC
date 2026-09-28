"""Independent, semi-analytic closed-flux Grad–Shafranov verification fixture."""
from dataclasses import dataclass
from functools import cached_property
import numpy as np
from numpy.polynomial import polynomial as poly
from scipy.constants import mu_0 as MU0
from scipy.optimize import brentq
from .mhd_benchmarks import radial_average
from .mhd_state import from_primitive


def _series(eigenvalue, axial_wavenumber, terms=48):
    """Regular solution F(s)=s² sum(a_n s^(2n)) of the radial GS ODE."""
    a = np.zeros(terms)
    a[0] = 1
    for n in range(1, terms):
        a[n] = (axial_wavenumber**2*a[n-1]-eigenvalue*(a[n-2] if n > 1 else 0))/(4*n*(n+1))
    return a


@dataclass(frozen=True)
class PeriodicEquilibrium:
    """Gas-filled, conducting-wall periodic equilibrium; NOT coil formation.

    psi = B0 R² F(r/R) cos(pi z/Z), p = p0 + lambda psi²/(2 mu0).
    F''-F'/s + (Lambda s²-k²) F=0, F(0)=F(1)=0, F/s² -> 1.
    lambda=Lambda/R⁴ and k=pi R/Z. The first radial eigenfunction has
    an internal O-point, reversed Bz and curved closed poloidal flux surfaces.
    """
    radius: float = 1.0
    half_length: float = 1.0
    field_scale: float = np.sqrt(MU0)
    background_pressure: float = 1.0
    density: float = 1.0

    def __post_init__(self):
        if any(not np.isfinite(v) or v <= 0 for v in (self.radius, self.half_length, self.field_scale, self.background_pressure, self.density)):
            raise ValueError("Equilibrium dimensions, field, pressure and density must be finite and positive")
        if not .5 <= self.half_length/self.radius <= 4:
            raise ValueError("This verification series is restricted to 0.5 <= Z/R <= 4")

    @property
    def k(self):
        return np.pi*self.radius/self.half_length

    @cached_property
    def eigenvalue(self):
        residual = lambda value: float(np.sum(_series(value, self.k)))
        # Small brackets find the first zero rather than a higher radial mode.
        previous = 0.0
        for upper in np.arange(5.0, 505.0, 5.0):
            if residual(previous)*residual(upper) < 0:
                return brentq(residual, previous, upper, xtol=1e-12)
            previous = upper
        raise RuntimeError("Could not bracket the fundamental radial GS eigenvalue")

    @cached_property
    def coefficients(self):
        return _series(self.eigenvalue, self.k)

    def radial_flux(self, s):
        s = np.asarray(s)
        return s*s*poly.polyval(s*s, self.coefficients)

    def fields(self, r, z):
        s, angle = np.asarray(r)/self.radius, np.pi*np.asarray(z)/self.half_length
        g = poly.polyval(s*s, self.coefficients)
        dg = poly.polyval(s*s, poly.polyder(self.coefficients))
        psi = self.field_scale*self.radius**2*s*s*g*np.cos(angle)
        br = self.field_scale*self.k*s*g*np.sin(angle)
        bz = 2*self.field_scale*(g+s*s*dg)*np.cos(angle)
        pressure = self.background_pressure+self.eigenvalue*psi**2/(2*MU0*self.radius**4)
        return psi, br, bz, pressure

    def state(self, mesh, gamma=5/3):
        if mesh.radius != self.radius or mesh.half_length != self.half_length:
            raise ValueError("Equilibrium and mesh dimensions must agree")
        psi, _, _, _ = self.fields(*mesh.corners)
        psi[[0, -1]] = 0  # Remove only series/root roundoff at exact boundaries.
        psi[:, -1] = psi[:, 0]
        radial = radial_average(mesh, lambda r: self.radial_flux(r/self.radius)**2)
        wave = np.pi/self.half_length
        axial = .5+np.diff(np.sin(2*wave*mesh.zf))/(4*wave*mesh.dz)
        pressure = self.background_pressure+self.eigenvalue*self.field_scale**2/(2*MU0)*radial[:, None]*axial[None, :]
        return from_primitive(mesh, psi, self.density, pressure, gamma=gamma)

    def ode_residual(self, s):
        """Continuous differential residual, independent of the MHD grid."""
        x = np.asarray(s)**2
        a = self.coefficients
        g, dg, ddg = poly.polyval(x, a), poly.polyval(x, poly.polyder(a)), poly.polyval(x, poly.polyder(a, 2))
        return x*(8*dg+4*x*ddg+(self.eigenvalue*x-self.k**2)*g)

    @property
    def reference_alfven_speed(self):
        return self.field_scale/np.sqrt(MU0*self.density)
