"""Continuum spherical induction mode including electron drift kinetic energy."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0, electron_mass as ME, elementary_charge as E
from .spherical_decay_reference import SphericalDecay


@dataclass(frozen=True)
class InertialSphere:
    radius: float = 1.
    diffusivity: float = .1
    electron_length: float = .2
    central_a: float = 1e-8
    mode: int = 1

    def __post_init__(self):
        if not np.isfinite(self.electron_length) or self.electron_length<0:
            raise ValueError("Finite nonnegative electron inertial length required")
        SphericalDecay(self.radius,self.diffusivity,self.central_a,self.mode)

    @property
    def wavenumber(self):
        return self.mode*np.pi/self.radius

    @property
    def energy_factor(self):
        return 1+(self.electron_length*self.wavenumber)**2

    @property
    def electron_density(self):
        return None if self.electron_length==0 else ME/(MU0*E**2*self.electron_length**2)

    @property
    def radial_reference(self):
        return SphericalDecay(self.radius,self.diffusivity/self.energy_factor,self.central_a,self.mode)

    def fields(self,r,z,time=0.):
        return self.radial_reference.fields(r,z,time)

    def heat_density(self,r,z,time):
        return self.energy_factor*self.radial_reference.heat_density(r,z,time)

    def integrals(self,time=0.):
        result=self.radial_reference.integrals(time)
        result["current_kinetic_energy_j"]=(self.energy_factor-1)*result["induced_energy_j"]
        result["total_response_energy_j"]=self.energy_factor*result["induced_energy_j"]
        for key in ("ohmic_power_w","ohmic_heat_j"):result[key]*=self.energy_factor
        return result
