"""Independent warm species extension of finite-sound-speed Hall eigenmodes."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import elementary_charge as E
from scipy.special import j0, j1
from .compressible_hall_reference import CompressibleHallWave


@dataclass(frozen=True)
class WarmHallWave(CompressibleHallWave):
    electron_fraction: float = 5/6

    def __post_init__(self):
        super().__post_init__()
        if not np.isfinite(self.electron_fraction) or not 0 <= self.electron_fraction < 1:
            raise ValueError('Electron fraction must be finite in [0,1)')

    @property
    def electron_pressure(self):
        return self.electron_fraction*self.pressure

    def fields(self, r, z, time=0., *, complex_output=False):
        out = super().fields(r, z, time, complex_output=True)
        r, z = np.broadcast_arrays(r, z)
        phase = np.exp(1j*(self.k*z-self.omega*time))
        sigma = self.eigenpair[1][3]
        amplitude = self.electron_fraction*self.density*self.sound_speed*sigma*phase
        gradient = np.stack((-self.q*amplitude*j1(self.q*r), np.zeros_like(amplitude),
                              1j*self.k*amplitude*j0(self.q*r)), axis=-1)
        out['electron_pressure'] = self.electron_fraction*out['pressure']
        out['ion_pressure'] = (1-self.electron_fraction)*out['pressure']
        out['electron_pressure_gradient'] = gradient
        out['electric'] -= gradient/(self.number_density*E)
        out['electron_velocity'] = out['velocity']-out['current']/(self.number_density*E)
        return out if complex_output else {k: v.real for k, v in out.items()}
