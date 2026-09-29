"""Nonlinear implicit midpoint Hall/resistive substep with explicit solver residuals."""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import block_diag, bmat, diags
from scipy.sparse.linalg import splu, gmres, LinearOperator
from scipy.constants import elementary_charge as E
from .hall_matrices import magnetic_hessian, transport_matrix
from .compatible_poynting import energy_rate


class HallConvergenceError(RuntimeError):
    """The nonlinear trial was not accepted; input arrays remain unchanged."""


@dataclass(frozen=True)
class HallMidpointResult:
    psi: np.ndarray
    g: np.ndarray
    total_energy_change: np.ndarray
    physical_heat_j: float
    iterations: int
    nonlinear_residual: float
    linear_residual: float
    electric_coefficients: np.ndarray


class HallMidpoint:
    def __init__(self, geometry, *, tolerance=1e-10, max_iterations=16):
        if not np.isfinite(tolerance) or not 0 < tolerance < 1 or type(max_iterations) is not int or max_iterations < 1:
            raise ValueError('Positive nonlinear tolerance below one and positive iteration budget required')
        self.geometry, self.tolerance, self.max_iterations = geometry, tolerance, max_iterations
        self.mass = block_diag(geometry.mass, format='csc')
        self.hessian = magnetic_hessian(geometry.mesh)
        self.stiffness = (geometry.faraday.T @ self.hessian @ geometry.faraday).tocsc()
        self.scale = 1/np.sqrt(self.mass.diagonal())
        self.diagonal = diags(self.scale)
        self.scaled_mass = (self.diagonal @ self.mass @ self.diagonal).tocsc()
        self._preconditioner = None
        self.factorizations = 0

    def _solve(self, matrix, rhs):
        if self._preconditioner is not None:
            preconditioner = LinearOperator(matrix.shape, matvec=self._preconditioner.solve, dtype=float)
            solved, info = gmres(matrix, rhs, M=preconditioner, restart=24, maxiter=3,
                                 rtol=max(4*np.finfo(float).eps, min(1e-12, .1*self.tolerance)), atol=0.)
            if info == 0:
                return solved
        self._preconditioner = splu(matrix)
        self.factorizations += 1
        return self._preconditioner.solve(rhs)

    def step(self, psi, g, number_density, dt, *, resistivity=0.):
        geometry = self.geometry
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Finite positive Hall step required')
        psi, g = np.asarray(psi, float), np.asarray(g, float)
        mesh = geometry.mesh
        if psi.shape != (mesh.nr+1, mesh.nz+1) or g.shape != (mesh.nr, mesh.nz) or not np.isfinite(psi).all() or not np.isfinite(g).all():
            raise ValueError('Finite corner psi and cell g required')
        flux_scale = max(float(np.max(abs(psi))), np.finfo(float).tiny)
        if np.max(abs(psi[0])) > 1e-13*flux_scale or np.ptp(psi[-1]) > 1e-13*flux_scale or np.max(abs(psi[:, 0]-psi[:, -1])) > 1e-13*flux_scale:
            raise ValueError('Regular axis, constant wall flux and periodic psi required')
        old = geometry.pack(psi, g).copy()
        if not np.isfinite(old).all():
            raise ValueError('Finite magnetic state required')
        eta = np.broadcast_to(np.asarray(resistivity, float), geometry.r.shape)
        density = None if number_density is None else np.broadcast_to(np.asarray(number_density, float), geometry.r.shape)
        magnetic = geometry.magnetic_samples(psi, g)
        transport_matrix(geometry, magnetic, density, eta)  # Validate before shortcuts.
        impedance = float(eta.max())
        if density is not None:
            impedance = max(impedance, float(np.max(np.linalg.norm(magnetic, axis=-1)/(density*E))))
        if impedance == 0:
            return HallMidpointResult(psi.copy(), g.copy(), np.zeros_like(g), 0., 0, 0., 0., np.zeros(self.mass.shape[0]))
        source = -geometry.faraday.T @ (self.hessian @ old)
        rhs = np.concatenate((self.scale*source, np.zeros_like(source)))
        upper = (.5*dt*impedance*self.diagonal @ self.stiffness @ self.diagonal).tocsc()
        candidate, linear_error = old.copy(), 0.
        for iteration in range(1, self.max_iterations+1):
            midpsi, midg = geometry.unpack(.5*(old+candidate))
            magnetic = geometry.magnetic_samples(midpsi, midg)
            transport = transport_matrix(geometry, magnetic, density, eta)
            lower = -(self.diagonal @ transport @ self.diagonal)/impedance
            matrix = bmat([[self.scaled_mass, upper], [lower, self.scaled_mass]], format='csc')
            solved = self._solve(matrix, rhs)
            residual = matrix @ solved-rhs
            linear_error = max(linear_error, float(np.linalg.norm(residual)/max(np.linalg.norm(rhs), np.finfo(float).tiny)))
            current = self.scale*solved[:len(source)]
            electric = impedance*self.scale*solved[len(source):]
            candidate = old+dt*(geometry.faraday @ electric)
            if not np.isfinite(candidate).all():
                raise HallConvergenceError('Nonfinite Hall midpoint trial; state not committed')
            midpsi, midg = geometry.unpack(.5*(old+candidate))
            actual_transport = transport_matrix(geometry, geometry.magnetic_samples(midpsi, midg), density, eta)
            electric_mass, electric_source = self.mass @ electric, actual_transport @ current
            scale = np.linalg.norm(self.scale*electric_mass)+np.linalg.norm(self.scale*electric_source)
            nonlinear_error = float(np.linalg.norm(self.scale*(electric_mass-electric_source))/max(scale, np.finfo(float).tiny))
            if nonlinear_error <= self.tolerance and linear_error <= max(self.tolerance, 1e-12):
                samples = geometry.sample_electric(current)
                heat = float(dt*np.sum(geometry.weight*eta*np.sum(samples*samples, axis=-1)))
                change = dt*energy_rate(geometry, midpsi, midg, electric)
                final_psi, final_g = geometry.unpack(candidate)
                return HallMidpointResult(final_psi, final_g, change, heat, iteration, nonlinear_error, linear_error, electric)
        raise HallConvergenceError(f'Hall midpoint failed after {self.max_iterations} iterations; relative residual {nonlinear_error:.6g}; state not committed')
