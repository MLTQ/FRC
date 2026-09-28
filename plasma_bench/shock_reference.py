"""Independent perpendicular-shock jumps and cylindrical steady-flow branches."""
from dataclasses import dataclass
from functools import cached_property
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.optimize import brentq


@dataclass(frozen=True)
class CylindricalShock:
    radius: float = .5
    upstream_density: float = 1.
    upstream_pressure: float = .1
    upstream_velocity: float = -6.
    upstream_field: float = float(np.sqrt(MU0))
    gamma: float = 5/3

    @cached_property
    def downstream(self):
        """Entropy-increasing compressive Rankine–Hugoniot root, not solver data."""
        rho, p, v, b, gamma = self.upstream_density, self.upstream_pressure, self.upstream_velocity, self.upstream_field, self.gamma
        if rho <= 0 or p <= 0 or gamma <= 1 or v >= 0 or self.radius <= 0:
            raise ValueError("Standing inward shock requires positive rho/p/radius and inward flow")
        va2, sound2 = b*b/(MU0*rho), gamma*p/rho
        if v*v <= sound2+va2:
            raise ValueError("Upstream perpendicular flow must be super-fast-magnetosonic")
        a = gamma/(gamma-1)
        enthalpy = a*p/rho+.5*v*v+va2
        polynomial = [(1-a/2)*va2, -enthalpy, a*(p/rho+v*v+.5*va2), (.5-a)*v*v]
        roots = np.roots(np.trim_zeros(polynomial, "f"))
        physical = []
        for root in roots:
            if abs(root.imag) > 1e-10 or root.real <= 1+1e-10:
                continue
            compression = float(root.real)
            pressure = p+rho*v*v*(1-1/compression)+b*b/(2*MU0)*(1-compression**2)
            if pressure > 0 and pressure/(rho*compression)**gamma > p/rho**gamma:
                physical.append((rho*compression, v/compression, pressure, b*compression))
        if len(physical) != 1:
            raise ValueError("Expected one admissible perpendicular shock root")
        return physical[0]

    @property
    def mass_flux(self):
        return self.radius*self.upstream_density*self.upstream_velocity

    @property
    def field_ratio(self):
        return self.upstream_field/self.upstream_density

    @property
    def specific_energy(self):
        return self.gamma/(self.gamma-1)*self.upstream_pressure/self.upstream_density+\
            .5*self.upstream_velocity**2+self.upstream_field**2/(MU0*self.upstream_density)

    def branch(self, radius, *, downstream):
        """Solve the cylindrical Bernoulli invariant on its sub/super-fast branch."""
        r = np.asarray(radius, float)
        if np.any(r <= 0) or not np.all(np.isfinite(r)):
            raise ValueError("Steady radial flow is defined only at finite r>0")
        rho0, _, p0, _ = self.downstream if downstream else (self.upstream_density, self.upstream_velocity, self.upstream_pressure, self.upstream_field)
        entropy = p0/rho0**self.gamma
        a = self.gamma/(self.gamma-1)
        magnetic = self.field_ratio**2/MU0
        density = np.empty_like(r)
        for index in np.ndindex(r.shape):
            m2 = (self.mass_flux/r[index])**2
            def residual(y):
                return a*entropy*np.exp((self.gamma-1)*y)+.5*m2*np.exp(-2*y)+magnetic*np.exp(y)-self.specific_energy
            def derivative(y):
                return self.gamma*entropy*np.exp((self.gamma-1)*y)-m2*np.exp(-2*y)+magnetic*np.exp(y)
            low, high = np.log(rho0)-25, np.log(rho0)+25
            critical = brentq(derivative, low, high, xtol=1e-14)
            if residual(critical) >= 0:
                raise ValueError(f"No steady branch beyond the cylindrical choking point at r={r[index]:g}")
            y = brentq(residual, critical if downstream else low, high if downstream else critical, xtol=1e-14)
            density[index] = np.exp(y)
        return density, self.mass_flux/(r*density), entropy*density**self.gamma, self.field_ratio*density

    def profile(self, radius):
        r = np.asarray(radius, float)
        result = [np.empty_like(r) for _ in range(4)]
        for down in (False, True):
            mask = (r < self.radius) if down else (r >= self.radius)
            if np.any(mask):
                values = self.branch(r[mask], downstream=down)
                for out, value in zip(result, values):
                    out[mask] = value
        return tuple(result)

    def flux_invariants(self, density, velocity, pressure, field):
        """Mass, momentum, energy and induction flux across a planar jump."""
        energy = pressure/(self.gamma-1)+.5*density*velocity**2+field**2/(2*MU0)
        return np.asarray([density*velocity, density*velocity**2+pressure+field**2/(2*MU0),
                           velocity*(energy+pressure+field**2/(2*MU0)), velocity*field])
