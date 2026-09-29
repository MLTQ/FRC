"""Exact decaying sphere for a prescribed variable carrier/relaxation profile."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0, electron_mass as ME, elementary_charge as E
from .spherical_decay_reference import SphericalDecay


@dataclass(frozen=True)
class ProfiledSphere:
    radius: float = 1.
    electron_length: float = .3
    profile_strength: float = 2.
    decay_rate: float = .7
    central_a: float = 1e-8

    def __post_init__(self):
        values=(self.radius,self.electron_length,self.profile_strength,self.decay_rate,self.central_a)
        if not all(np.isfinite(v) for v in values) or min(self.radius,self.electron_length,self.decay_rate)<=0 or self.profile_strength<0:
            raise ValueError("Finite positive scales and nonnegative profile strength required")

    def de_squared(self,r,z):
        return self.electron_length**2*(1+self.profile_strength*(r*r+z*z)/self.radius**2)

    def density(self,r,z):
        return ME/(MU0*E**2*self.de_squared(r,z))

    def resistivity(self,r,z):
        return MU0*self.decay_rate*((self.radius/np.pi)**2+self.de_squared(r,z))

    @property
    def radial_reference(self):
        return SphericalDecay(self.radius,self.decay_rate*(self.radius/np.pi)**2,self.central_a)

    def fields(self,r,z,time=0.):
        return self.radial_reference.fields(r,z,time)

    def heat_density(self,r,z,time):
        if not np.isfinite(time) or time<0:raise ValueError("Finite nonnegative reference time required")
        j=self.fields(r,z,0.)['Jtheta']
        return self.resistivity(r,z)*j*j*(-np.expm1(-2*self.decay_rate*time))/(2*self.decay_rate)
