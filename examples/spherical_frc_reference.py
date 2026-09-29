"""Analytic ideal-MHD spherical FRC with a current-free infinite exterior.

This prescribed equilibrium is a future free-boundary verification reference,
not a formation calculation, density closure, or stability/confinement claim.
"""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0


@dataclass(frozen=True)
class SphericalFRCReference:
    radius: float = 1.
    bias: float = np.sqrt(MU0)

    def __post_init__(self):
        if not np.isfinite(self.radius) or self.radius <= 0 or not np.isfinite(self.bias) or self.bias == 0:
            raise ValueError("Positive finite radius and finite nonzero bias required")

    def fields(self,r,z):
        """Return SI psi, Br, Bz, pressure and Jtheta; no exterior fluid."""
        r,z=np.broadcast_arrays(np.asarray(r,dtype=float),np.asarray(z,dtype=float))
        if np.any(r<0) or not np.all(np.isfinite(r+z)):
            raise ValueError("Finite coordinates and r >= 0 required")
        a,b=self.radius,self.bias
        c=-3*b/4;dipole=-b*a**3/2
        s2=r*r+z*z
        inside=s2<=a*a
        safe_s=np.sqrt(np.maximum(s2,a*a))
        psi=np.where(inside,c*r*r*(1-s2/(a*a)),.5*b*r*r+dipole*r*r/safe_s**3)
        br=np.where(inside,2*c*r*z/(a*a),3*dipole*r*z/safe_s**5)
        bz=np.where(inside,2*c*(1-(2*r*r+z*z)/(a*a)),b+dipole*(2*z*z-r*r)/safe_s**5)
        pressure=np.where(inside,10*c*psi/(MU0*a*a),0.)
        current=np.where(inside,10*c*r/(MU0*a*a),0.)
        return psi,br,bz,pressure,current

    def integrals(self,gamma=5/3):
        """Independent continuum integrals; uniform-bias self energy excluded."""
        if not np.isfinite(gamma) or gamma<=1:
            raise ValueError("gamma must be finite and greater than one")
        a,b=self.radius,self.bias
        scale=np.pi*b*b*a**3/MU0
        return {"separatrix_volume_m3":4*np.pi*a**3/3,
                "magnetic_axis_r_m":a/np.sqrt(2),"axis_Bz_over_bias":-1.5,
                "closed_poloidal_flux_wb":3*np.pi*abs(b)*a*a/8,
                "pressure_integral_j":6*scale/7,
                "internal_energy_j":6*scale/(7*(gamma-1)),
                "inside_total_magnetic_energy_j":3*scale/7,
                "induced_field_energy_j":10*scale/7,
                "exterior_induced_field_energy_j":scale/3,
                "toroidal_current_a":-5*b*a/MU0}
