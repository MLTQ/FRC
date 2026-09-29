"""Independent sinusoidal diffusion solution for a stationary infinite cylinder."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.special import ive


@dataclass(frozen=True)
class CylindricalSkin:
    radius: float = 1.
    diffusivity: float = .1
    frequency: float = 10.
    field_amplitude: float = 1e-4*np.sqrt(MU0)

    def __post_init__(self):
        if any(not np.isfinite(v) or v <= 0 for v in
               (self.radius, self.diffusivity, self.frequency, self.field_amplitude)):
            raise ValueError("Skin-reference parameters must be finite and positive")

    @property
    def omega(self):
        return 2*np.pi*self.frequency

    @property
    def skin_depth(self):
        return np.sqrt(2*self.diffusivity/self.omega)

    def phasors(self, radius):
        """Return psi, Bz and Jtheta, with real(field * exp(i omega t))."""
        r = np.asarray(radius, float)
        if np.any(~np.isfinite(r)) or np.any(r < 0):
            raise ValueError("Evaluation radii must be finite and nonnegative")
        k = (1+1j)/self.skin_depth
        inside = np.minimum(r, self.radius)
        scale = np.exp(k.real*(inside-self.radius))/ive(0, k*self.radius)
        b = self.field_amplitude*ive(0, k*inside)*scale
        j = -self.field_amplitude*k*ive(1, k*inside)*scale/MU0
        psi = self.field_amplitude*inside*ive(1, k*inside)*scale/k
        psi = psi+.5*self.field_amplitude*np.maximum(r*r-self.radius**2, 0)
        return psi, np.where(r <= self.radius, b, self.field_amplitude), np.where(r <= self.radius, j, 0.)

    def fields(self, radius, time):
        return tuple(np.real(value*np.exp(1j*self.omega*time)) for value in self.phasors(radius))

    def mean_ohmic_power_per_length(self):
        """Independent complex Poynting flux at r=R, positive inward [W/m]."""
        _, b, j = self.phasors(self.radius)
        electric = self.diffusivity*MU0*j
        return float(-np.pi*self.radius*np.real(electric*np.conj(b))/MU0)
