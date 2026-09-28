"""Stationary resistive induction coupled to an exactly nonconducting vacuum."""
from dataclasses import dataclass
from typing import Protocol
import numpy as np
from scipy.constants import mu_0 as MU0
from scipy.sparse import diags
from scipy.sparse.linalg import splu
from .geometry import Grid
from .numerics import stiffness, expand


class PlasmaSolver(Protocol):
    def step(self, external_psi, dt: float): ...


@dataclass
class DiffusionStep:
    induced_psi: np.ndarray
    e_theta: np.ndarray
    j_theta: np.ndarray
    induced_energy: float
    ohmic_heat: float
    source_work: float
    numerical_loss: float
    balance_residual: float


class StationaryDiffusion:
    def __init__(self, grid: Grid, conductivity, initial_external_psi, initial_induced=None):
        self.grid = grid
        self.sigma = np.broadcast_to(np.asarray(conductivity, dtype=float), grid.shape).copy()
        if not np.all(np.isfinite(self.sigma)) or np.any(self.sigma < 0):
            raise ValueError("Conductivity must be finite and nonnegative")
        self.k = stiffness(grid)
        self.r = grid.r[1:-1, None]
        self.mass = (MU0*self.sigma[1:-1, 1:-1]*grid.dr*grid.dz/self.r).ravel()
        self.external = self._check_flux(initial_external_psi).copy()
        self.u = np.zeros_like(self.mass) if initial_induced is None else self._check_flux(initial_induced)[1:-1, 1:-1].ravel().copy()
        if initial_induced is not None:
            initial = np.asarray(initial_induced)
            if np.any(initial[[0, -1], :]) or np.any(initial[:, [0, -1]]):
                raise ValueError("Initial induced flux must vanish at boundaries")
            ku = self.k@self.u
            if np.linalg.norm(ku[self.mass == 0]) > 1e-10*max(np.linalg.norm(ku), 1e-30):
                raise ValueError("Initial induced field violates the vacuum constraint")
        self._factor = None
        self._dt = None

    def _check_flux(self, flux):
        p = np.asarray(flux, dtype=float)
        if p.shape != self.grid.shape or not np.all(np.isfinite(p)):
            raise ValueError("External/initial flux must be finite with grid shape")
        if np.any(np.abs(p[0]) > 1e-25):
            raise ValueError("Flux must use the zero-axis gauge")
        return p

    @property
    def induced_energy(self):
        return float(np.pi/MU0*(self.u@(self.k@self.u)))

    def step(self, external_psi, dt: float):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("dt must be finite and positive")
        external = self._check_flux(external_psi)
        if self._dt is None or not np.isclose(dt, self._dt, rtol=1e-12, atol=0):
            self._factor = splu(diags(self.mass, format="csc")+dt*self.k)
            self._dt = dt
        dp = external-self.external
        old = self.u
        before = self.induced_energy
        self.u = self._factor.solve(self.mass*(old-dp[1:-1, 1:-1].ravel()))
        du = self.u-old
        total_delta = dp+expand(self.grid, du)
        r, _ = self.grid.mesh
        electric = np.divide(-total_delta/dt, r, out=np.zeros_like(r), where=r > 0)
        j = self.sigma*electric
        cell_area = self.grid.dr*self.grid.dz
        active_j = j[1:-1, 1:-1]
        # Equivalent sigma*E^2 expression handles sigma=0 without division.
        heat = float(dt*2*np.pi*cell_area*np.sum(self.r*self.sigma[1:-1, 1:-1]*electric[1:-1, 1:-1]**2))
        work = float(-2*np.pi*cell_area*np.sum(active_j*dp[1:-1, 1:-1]))
        loss = float(np.pi/MU0*(du@(self.k@du)))
        energy = self.induced_energy
        residual = energy-before+heat+loss-work
        self.external = external.copy()
        return DiffusionStep(expand(self.grid, self.u), electric, j, energy, heat, work, loss, residual)


class ResistiveMHD:
    """Formation-stage gate until dynamic plasma/vacuum coupling is verified."""
    def __init__(self, *args, **kwargs):
        raise NotImplementedError("Stage 2 formation requires verified dynamic plasma/vacuum coupling. The developmental conservative kernel is plasma_bench.mhd.AxisymmetricMHD; its conducting-wall/periodic verification boundaries are not the quartz-tube experiment. See docs/MHD_NUMERICS.md.")
