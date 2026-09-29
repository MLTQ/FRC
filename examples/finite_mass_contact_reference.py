"""Exact unforced finite-mass two-temperature contact, independent of solvers."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import atomic_mass,electron_mass,elementary_charge as E,mu_0 as MU0


@dataclass(frozen=True)
class FiniteMassContact:
    radius:float=.035
    length:float=.24
    number_density:float=1e20
    ion_mass:float=39.948*atomic_mass
    electron_mass:float=electron_mass
    gamma_e:float=5/3
    gamma_i:float=1.4
    electron_ev:float=5.
    ion_ev:float=1.
    density_amplitude:float=.25
    pressure_amplitude:float=.12
    velocity:float=2400.
    bias:float=.005

    @property
    def period(self):return self.length/abs(self.velocity)

    @property
    def total_pressure(self):return self.number_density*E*(self.electron_ev+self.ion_ev)

    def fields(self,z,time=0.,*,rate=False):
        k=2*np.pi/self.length;phase=k*(z-self.velocity*time);pe0=self.number_density*E*self.electron_ev
        n=self.number_density*(1+self.density_amplitude*np.cos(phase));pe=pe0*(1+self.pressure_amplitude*np.sin(phase+.23))
        pi=self.total_pressure-pe;rho=(self.ion_mass+self.electron_mass)*n
        gradient=pe0*self.pressure_amplitude*k*np.cos(phase+.23)
        if rate:
            n=self.number_density*self.density_amplitude*k*self.velocity*np.sin(phase)
            rho=(self.ion_mass+self.electron_mass)*n;pe=-self.velocity*gradient;pi=-pe
        energy=.5*rho*self.velocity**2+pe/(self.gamma_e-1)+pi/(self.gamma_i-1)
        if not rate:energy+=self.bias**2/(2*MU0)
        result=dict(number_density=n,density=rho,momentum=rho*self.velocity,electron_pressure=pe,ion_pressure=pi,energy=energy)
        if not rate:result['electric_z']=-gradient/(E*n)
        return result
