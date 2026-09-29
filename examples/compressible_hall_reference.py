"""Independent finite-sound-speed, cold-electron cylindrical Hall wave oracle."""
from dataclasses import dataclass
from functools import cached_property
import numpy as np
from scipy.constants import elementary_charge as E, mu_0 as MU0
from scipy.linalg import eigh
from scipy.special import j0, j1, jn_zeros


def generator(q, k, sound_speed, alfven_speed, hall_diffusivity):
    """Continuum generator in equal-energy coordinates (u,w,v,sigma,bp,bt).

    sigma=p1/(rho cs), bp=K a/sqrt(mu0 rho), bt=Btheta/sqrt(mu0 rho).
    Signed vA and DH retain the sign of the axial bias. q=0 is a separate
    Cartesian parallel limit, not a nontrivial cylindrical J1 eigenfunction.
    """
    if not np.isfinite([q, k, sound_speed, alfven_speed, hall_diffusivity]).all() or q < 0 or sound_speed <= 0:
        raise ValueError('Finite coefficients, nonnegative radial wavenumber and positive sound speed required')
    total = np.hypot(q, k)
    a = np.zeros((6, 6), complex)
    a[0, 3], a[3, 0] = q*sound_speed, -q*sound_speed
    a[0, 4], a[4, 0] = total*alfven_speed, -total*alfven_speed
    a[1, 5] = a[5, 1] = 1j*k*alfven_speed
    a[2, 3] = a[3, 2] = -1j*k*sound_speed
    a[4, 5] = a[5, 4] = -1j*k*total*hall_diffusivity
    return a


def dispersion_terms(omega, q, k, sound_speed, alfven_speed, hall_diffusivity):
    """Independent eliminated sixth-degree dispersion: return LHS and RHS."""
    x, total2 = np.asarray(omega)**2, q*q+k*k
    lhs = (x-k*k*alfven_speed**2)*(x*x-total2*(sound_speed**2+alfven_speed**2)*x
                                +k*k*total2*alfven_speed**2*sound_speed**2)
    rhs = x*hall_diffusivity**2*k*k*total2*(x-total2*sound_speed**2)
    return lhs, rhs


@dataclass(frozen=True)
class CompressibleHallWave:
    radius: float
    length: float
    number_density: float
    ion_mass: float
    bias: float
    pressure: float
    gamma: float = 5/3
    radial_mode: int = 1
    axial_mode: int = 1
    branch: int = 2  # ascending positive frequency; not a fixed helicity label
    sense: int = 1   # sign of omega in exp(i k z-i omega t)
    amplitude: float = 1e-4  # equal-energy coordinate norm, m/s
    hall: bool = True

    def __post_init__(self):
        if not np.isfinite([self.radius, self.length, self.number_density, self.ion_mass,
                            self.bias, self.pressure, self.gamma, self.amplitude]).all():
            raise ValueError('Finite wave inputs required')
        if min(self.radius, self.length, self.number_density, self.ion_mass, self.pressure) <= 0 or self.gamma <= 1 or self.bias == 0 or self.amplitude <= 0:
            raise ValueError('Positive dimensions/density/mass/pressure/amplitude, gamma > 1 and nonzero bias required')
        if type(self.radial_mode) is not int or self.radial_mode < 1 or type(self.axial_mode) is not int or self.axial_mode == 0:
            raise ValueError('Positive integer radial mode and nonzero integer axial mode required')
        if type(self.branch) is not int or self.branch not in (0, 1, 2) or self.sense not in (-1, 1) or type(self.hall) is not bool:
            raise ValueError('Branch 0, 1 or 2, sense +/-1 and boolean Hall selection required')

    @cached_property
    def q(self):
        return jn_zeros(1, self.radial_mode)[-1]/self.radius

    @property
    def k(self):
        return 2*np.pi*self.axial_mode/self.length

    @property
    def total_k(self):
        return np.hypot(self.q, self.k)

    @property
    def density(self):
        return self.number_density*self.ion_mass

    @property
    def sound_speed(self):
        return np.sqrt(self.gamma*self.pressure/self.density)

    @property
    def alfven_speed(self):
        return self.bias/np.sqrt(MU0*self.density)

    @property
    def hall_diffusivity(self):
        return self.bias/(MU0*self.number_density*E) if self.hall else 0.

    @cached_property
    def eigenpair(self):
        matrix = generator(self.q, self.k, self.sound_speed, self.alfven_speed, self.hall_diffusivity)
        # i A is Hermitian, with eigenvalue omega for A y=-i omega y.
        values, vectors = eigh(1j*matrix)
        index = 3+self.branch if self.sense == 1 else 2-self.branch
        vector = vectors[:, index]
        anchor = np.argmax(abs(vector))
        vector = self.amplitude*vector*np.exp(-1j*np.angle(vector[anchor]))
        vector.setflags(write=False)
        return float(values[index]), vector

    @property
    def omega(self):
        return self.eigenpair[0]

    def fields(self, r, z, time=0., *, complex_output=False):
        """Perturbations only: vectors ordered (r,theta,z), plus psi perturbation."""
        r, z = np.broadcast_arrays(r, z)
        phase = np.exp(1j*(self.k*z-self.omega*time))
        u, w, v, sigma, bp, bt = self.eigenpair[1]
        a, t = np.sqrt(MU0*self.density)*bp/self.total_k, np.sqrt(MU0*self.density)*bt
        f0, f1 = j0(self.q*r), j1(self.q*r)
        pressure = self.density*self.sound_speed*sigma*f0*phase
        magnetic = np.stack((-1j*self.k*a*f1*phase, t*f1*phase, self.q*a*f0*phase), axis=-1)
        velocity = np.stack((u*f1*phase, w*f1*phase, v*f0*phase), axis=-1)
        current = np.stack((-1j*self.k*t*f1*phase, self.total_k**2*a*f1*phase,
                            self.q*t*f0*phase), axis=-1)/MU0
        bias = np.array([0., 0., self.bias])
        electric = -np.cross(velocity, bias)
        if self.hall:
            electric += np.cross(current, bias)/(self.number_density*E)
        result = dict(magnetic=magnetic, velocity=velocity, current=current, electric=electric,
                      pressure=pressure, density=pressure/self.sound_speed**2, psi=r*a*f1*phase)
        return result if complex_output else {key: value.real for key, value in result.items()}

    def perturbation_energy(self):
        """Continuum quadratic wave energy, real fields over one periodic cylinder."""
        radial_norm = .5*self.radius**2*j0(self.q*self.radius)**2
        return .5*np.pi*self.length*radial_norm*self.density*self.amplitude**2
