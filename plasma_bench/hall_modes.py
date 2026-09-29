"""Work-paired poloidal/toroidal linear Hall dynamics and midpoint integration."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from .hall_parameters import HallParameters


@dataclass(frozen=True)
class HallUpdate:
    state: np.ndarray
    midpoint: np.ndarray
    heat: float
    balance_residual: float
    equation_relative_residual: float


class LinearHallModes:
    """Mode coefficients [Bp, Bt, sqrt(mu rho) Vp, sqrt(mu rho) Vt].

    Each poloidal basis vector is normalized to its toroidal partner. Curl
    exchanges those partners with coefficient sqrt(k_perp**2+k_z**2).
    Energy uses the caller's physical volume normalization / (2 mu0).
    """
    def __init__(self, kz, radial_eigenvalues, parameters: HallParameters):
        self.parameters = parameters
        kz = np.asarray(kz, dtype=float)
        radial = np.asarray(radial_eigenvalues, dtype=float)
        if kz.ndim != 1 or radial.ndim != 1 or not kz.size or not radial.size:
            raise ValueError('Wave numbers and radial eigenvalues must be nonempty vectors')
        if not np.isfinite(kz).all() or not np.isfinite(radial).all() or np.any(radial < 0):
            raise ValueError('Wave numbers finite; radial eigenvalues nonnegative')
        self.kz = kz.copy()
        self.k_squared = radial[:, None]+kz[None, :]**2
        self.shape = self.k_squared.shape+(4,)
        self.mass = np.ones(self.shape)
        self.mass[..., :2] += parameters.current_length_squared*self.k_squared[..., None]
        self.operator = np.zeros(self.shape+(4,), complex)
        damping = parameters.resistivity/MU0*self.k_squared
        hall = parameters.hall_coefficient*parameters.bias/MU0*kz*np.sqrt(self.k_squared)
        coupling = 1j*kz*parameters.alfven
        for j in (0, 1):
            self.operator[..., j, j] = -damping
            self.operator[..., j, 1-j] = -1j*hall
            self.operator[..., j, j+2] = coupling
            self.operator[..., j+2, j] = coupling
        self._cached_dt = None

    def validate(self, state):
        state = np.asarray(state, dtype=complex)
        if state.shape != self.shape or not np.isfinite(state).all():
            raise ValueError(f'Finite mode state with shape {self.shape} required')
        return state

    def energies(self, state, volume_factor=1.):
        state = self.validate(state)
        if not np.isfinite(volume_factor) or volume_factor <= 0:
            raise ValueError('Positive finite volume factor required')
        field = np.sum(abs(state[..., :2])**2, axis=-1)
        scale = volume_factor/(2*MU0)
        return {'magnetic': float(scale*np.sum(field)),
                'bulk_kinetic': float(scale*np.sum(abs(state[..., 2:])**2)),
                'current_kinetic': float(scale*self.parameters.current_length_squared
                                         *np.sum(self.k_squared*field))}

    def step(self, state, dt, volume_factor=1.):
        state = self.validate(state)
        if not np.isfinite(dt) or dt <= 0 or not np.isfinite(volume_factor) or volume_factor <= 0:
            raise ValueError('Positive finite time step and volume factor required')
        if self._cached_dt != dt:
            identity = np.eye(4)*self.mass[..., None]
            # Only the current step size is cached: no unbounded cache during adaptation.
            self._left = identity-.5*dt*self.operator
            self._cached_dt = dt
        rhs = dt*np.einsum('...ij,...j->...i', self.operator, state)
        increment = np.linalg.solve(self._left, rhs[..., None])[..., 0]
        final = state+increment
        if not np.isfinite(final).all():
            raise FloatingPointError('Nonfinite linear Hall trial')
        middle = state+.5*increment
        power = self.parameters.resistivity/MU0**2*volume_factor*np.sum(
            self.k_squared*np.sum(abs(middle[..., :2])**2, axis=-1))
        before = sum(self.energies(state, volume_factor).values())
        after = sum(self.energies(final, volume_factor).values())
        residual = np.einsum('...ij,...j->...i', self._left, increment)-rhs
        denom = np.linalg.norm(rhs)+np.linalg.norm(self._left)*np.linalg.norm(increment)
        return HallUpdate(final, middle, float(dt*power), float(after-before+dt*power),
                          float(np.linalg.norm(residual)/max(denom, np.finfo(float).tiny)))
