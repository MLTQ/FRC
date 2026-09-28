"""Analytic loop fields and reusable prescribed vacuum field bases."""
from dataclasses import dataclass
import numpy as np
from scipy.special import ellipk, ellipe
from scipy.constants import mu_0 as MU0
from .coils import Coil
from .geometry import Grid


@dataclass
class Field:
    psi: np.ndarray
    br: np.ndarray
    bz: np.ndarray
    dpsi_dt: np.ndarray
    e_theta: np.ndarray

    @property
    def magnitude(self):
        return np.hypot(self.br, self.bz)


def loop_field(r, z, radius: float, z0: float = 0.0):
    """Return (psi, Br, Bz) per ampere-turn, broadcasting r and z."""
    r, h = np.broadcast_arrays(np.asarray(r, dtype=float), np.asarray(z, dtype=float)-z0)
    if not np.isfinite(radius) or radius <= 0 or not np.isfinite(z0) or np.any(r < 0) or not np.all(np.isfinite(r+h)):
        raise ValueError("Finite coordinates, r >= 0 and positive coil radius required")
    a = radius
    d = (a+r)**2 + h*h
    q = (a-r)**2 + h*h
    if np.any(q/d < 1e-14):
        raise ValueError("Evaluation on/too near an ideal filament is singular")
    m = 4*a*r/d
    k, e = ellipk(m), ellipe(m)
    # Cancellation-free expansion of (1-m/2) K-E through m^6.
    series = np.pi/2*m*m*(1/16 + m*(3/64 + m*(75/2048 + m*(245/8192 + m*6615/262144))))
    factor = np.where(m < 1e-3, series, (1-m/2)*k-e)
    psi = MU0*np.sqrt(d)/(2*np.pi)*factor
    safe_r = np.where(r == 0, 1, r)
    br = MU0*h/(2*np.pi*safe_r*np.sqrt(d))*(-k+(a*a+r*r+h*h)/q*e)
    bz = MU0/(2*np.pi*np.sqrt(d))*(k+(a*a-r*r-h*h)/q*e)
    # Axis series through r^2 in Bz, r^3 in Br; relative errors O((r/L)^4).
    s = a*a+h*h
    c = MU0*a*a/2
    b0 = c/s**1.5
    b1 = -3*c*h/s**2.5
    b2 = c*(12*h*h-3*a*a)/s**3.5
    b3 = c*h*(45*a*a-60*h*h)/s**4.5
    near = r/np.sqrt(s) < 1e-3
    br = np.where(near, -r*b1/2+r**3*b3/16, br)
    bz = np.where(near, b0-r*r*b2/4, bz)
    return psi, br, bz


class VacuumField:
    def __init__(self, grid: Grid, coils: tuple[Coil, ...], bias: float = 0.0):
        if not np.isfinite(bias):
            raise ValueError("Bias field must be finite")
        self.grid, self.coils, self.bias = grid, tuple(coils), bias
        self.r, self.z = grid.mesh
        self.bases = [loop_field(self.r, self.z, c.radius, c.z) for c in coils]

    def at(self, t: float):
        if not np.isfinite(t):
            raise ValueError("Time must be finite")
        psi = self.bias*self.r**2/2
        br = np.zeros(self.grid.shape)
        bz = np.full(self.grid.shape, self.bias, dtype=float)
        rate = np.zeros(self.grid.shape)
        for coil, (p, b_r, b_z) in zip(self.coils, self.bases):
            current, derivative = coil.current(t), coil.rate(t)
            if not np.all(np.isfinite([current, derivative])):
                raise ValueError("Waveform returned nonfinite source data")
            psi = psi + current*p
            br += current*b_r
            bz += current*b_z
            rate += derivative*p
        electric = np.divide(-rate, self.r, out=np.zeros_like(rate), where=self.r > 0)
        return Field(psi, br, bz, rate, electric)
