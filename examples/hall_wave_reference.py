"""Independent continuum cylindrical/parallel two-fluid wave references."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E, electron_mass as ME, mu_0 as MU0
from scipy.special import j0, j1, jn_zeros


def frequency(n, mi, bias, eta, k, q, helicity, model='two_fluid'):
    """Positive-real-frequency root for exp(i kz-i omega t), k != 0."""
    if helicity not in (-1, 1) or k == 0 or bias == 0:
        raise ValueError('A propagating reference needs helicity +/-1, nonzero k and bias')
    me = ME if model == 'two_fluid' else 0.
    rho = n*(mi+me)
    length2 = mi*me/((mi+me)*MU0*n*E**2)
    alpha = (mi-me)/((mi+me)*n*E) if model != 'mhd' else 0.
    total2 = q*q+k*k
    a = 1+length2*total2
    coupling2 = (k*bias)**2/(MU0*rho)
    c = helicity*alpha*bias*k*np.sqrt(total2)/MU0-1j*eta/MU0*total2
    root = np.sqrt(c*c+4*a*coupling2)
    # Avoid subtracting nearly equal real parts on the ion branch.
    return (c+root)/(2*a) if c.real >= 0 else 2*coupling2/(root-c)


@dataclass(frozen=True)
class CylindricalWave:
    radius: float
    length: float
    number_density: float
    ion_mass: float
    bias: float
    resistivity: float
    helicity: int
    model: str = 'two_fluid'
    radial_mode: int = 1
    axial_mode: int = 1
    amplitude: float = 1e-8  # vector-potential amplitude, T m

    @property
    def q(self):
        return jn_zeros(1, self.radial_mode)[-1]/self.radius

    @property
    def k(self):
        return 2*np.pi*self.axial_mode/self.length

    @property
    def total_k(self):
        return np.hypot(self.q, self.k)

    @property
    def omega(self):
        return frequency(self.number_density, self.ion_mass, self.bias, self.resistivity,
                         self.k, self.q, self.helicity, self.model)

    @property
    def density(self):
        return self.number_density*(self.ion_mass+(ME if self.model == 'two_fluid' else 0.))

    @property
    def velocity_ratio(self):
        return -self.k*self.bias/(MU0*self.density*self.omega)

    def potentials(self, r, z, time=0.):
        a = self.amplitude*j1(self.q*np.asarray(r))*np.exp(1j*(self.k*np.asarray(z)-self.omega*time))
        return a, self.helicity*self.total_k*a, self.velocity_ratio*a, self.velocity_ratio*self.helicity*self.total_k*a

    def fields(self, r, z, time=0., quantity='magnetic', complex_output=False):
        r, z = np.broadcast_arrays(r, z)
        phase = self.amplitude*np.exp(1j*(self.k*z-self.omega*time))
        field = np.stack((-1j*self.k*j1(self.q*r)*phase,
                          self.helicity*self.total_k*j1(self.q*r)*phase,
                          self.q*j0(self.q*r)*phase), axis=-1)
        if quantity == 'current':
            field *= self.helicity*self.total_k/MU0
        elif quantity == 'velocity':
            field *= self.velocity_ratio
        elif quantity != 'magnetic':
            raise ValueError('Unknown reference field')
        return field if complex_output else field.real

    def heat_density(self, r, z, time):
        """Analytic integral of eta |Re(J0 exp(-i omega t))|² at each point."""
        current = self.fields(r, z, quantity='current', complex_output=True)
        omega = self.omega
        def integral(rate):
            return time if rate == 0 else np.expm1(rate*time)/rate
        constant = np.sum(abs(current)**2, axis=-1)*integral(2*omega.imag)
        rotating = np.sum(current**2, axis=-1)*integral(-2j*omega)
        return self.resistivity*.5*(constant+rotating.real)
