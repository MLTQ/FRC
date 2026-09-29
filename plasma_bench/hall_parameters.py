"""SI coefficients for the linear, uniform, incompressible two-species limit."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, mu_0 as MU0


@dataclass(frozen=True)
class HallParameters:
    number_density: float
    ion_mass: float
    bias: float
    resistivity: float = 0.
    model: str = 'two_fluid'

    def __post_init__(self):
        if not all(np.isfinite(v) and v > 0 for v in (self.number_density, self.ion_mass)):
            raise ValueError('Uniform carrier density and ion mass must be positive and finite')
        if not np.isfinite(self.bias) or not np.isfinite(self.resistivity) or self.resistivity < 0:
            raise ValueError('Bias must be finite and resistivity finite and nonnegative')
        if self.model not in ('two_fluid', 'hall', 'mhd'):
            raise ValueError('model must explicitly select two_fluid, hall, or mhd')

    @property
    def electron_mass(self):
        return ME if self.model == 'two_fluid' else 0.

    @property
    def density(self):
        return self.number_density*(self.ion_mass+self.electron_mass)

    @property
    def alfven(self):
        """Signed background-field coupling, not the nonnegative wave speed."""
        return self.bias/np.sqrt(MU0*self.density)

    @property
    def hall_coefficient(self):
        if self.model == 'mhd':
            return 0.
        return ((self.ion_mass-self.electron_mass)/(self.ion_mass+self.electron_mass)
                /(E*self.number_density))

    @property
    def current_length_squared(self):
        reduced_mass = self.ion_mass*self.electron_mass/(self.ion_mass+self.electron_mass)
        return reduced_mass/(MU0*self.number_density*E**2)
