"""Independent continuum targets for the bounded external-engine comparison."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import c, elementary_charge, mu_0, proton_mass


@dataclass(frozen=True)
class ParallelWave:
    """Cold proton, massless-electron circular wave; SI, periodic z."""
    density: float = 1e18
    bias: float = 1e-3
    amplitude: float = 1e-7
    k_di: float = 1.
    helicity: int = 1

    @property
    def ion_length(self):
        return np.sqrt(proton_mass/(mu_0*self.density*elementary_charge**2))

    @property
    def k(self):
        return self.k_di/self.ion_length

    @property
    def omega_i(self):
        return elementary_charge*self.bias/proton_mass

    @property
    def omega(self):
        # omega² - h k² vA di omega - k² vA² = 0.
        root = np.sqrt(self.k_di**4+4*self.k_di**2)
        scaled = ((self.k_di**2+root)/2 if self.helicity == 1
                  else 2*self.k_di**2/(root+self.k_di**2))
        return self.omega_i*scaled

    @property
    def length(self):
        return 2*np.pi/self.k

    @property
    def velocity_ratio(self):
        return -self.k*self.bias/(mu_0*proton_mass*self.density*self.omega)

    def magnetic(self, z, time):
        phase = self.k*np.asarray(z)-self.omega*time
        return np.array([self.amplitude*np.cos(phase),
                         -self.helicity*self.amplitude*np.sin(phase)])


def relative_error(actual, exact):
    return float(np.linalg.norm(actual-exact)/np.linalg.norm(exact))


def mode_amplitude(values, phase):
    """Complex cosine coefficient on the actual field-component coordinates."""
    return complex(2*np.mean(np.asarray(values)*np.exp(-1j*phase)))


def vacuum_field(z, time, length=1., amplitude=1e-7, electric=False):
    """Ex=c By, a right-moving vacuum Maxwell wave with no particles."""
    result = amplitude*np.cos(2*np.pi*(np.asarray(z)-c*time)/length)
    return c*result if electric else result
