"""Joint magnetic/current-inertia/species midpoint with atomic returned trials."""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import bmat
from scipy.constants import elementary_charge as E
from .hall_midpoint import HallMidpoint, HallConvergenceError
from .hall_matrices import transport_matrix
from .thermal_means import logarithmic_mean
from .current_thermal_flux import flux


@dataclass(frozen=True)
class CurrentThermalStep:
    psi: np.ndarray
    g: np.ndarray
    electron_entropy: np.ndarray
    ion_entropy: np.ndarray
    total_energy_change: np.ndarray
    electric_coefficients: np.ndarray
    radial_energy_flux: np.ndarray
    axial_energy_flux: np.ndarray
    iterations: int
    nonlinear_residual: float
    linear_residual: float
    stored_state_residual: float


class CurrentThermalMidpoint(HallMidpoint):
    """Reuse only compatible matrix scaling and the checked sparse linear solve."""
    def __init__(self, model, *, tolerance=1e-10, max_iterations=80):
        super().__init__(model.geometry, tolerance=tolerance, max_iterations=max_iterations)
        self.model = model

    def step(self, psi, g, electron_entropy, ion_entropy, dt):
        model, geo, mesh = self.model, self.geometry, self.geometry.mesh
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Finite positive current step required')
        psi, g = np.asarray(psi, float), np.asarray(g, float)
        if psi.shape != (mesh.nr+1, mesh.nz+1) or g.shape != (mesh.nr, mesh.nz) or not np.isfinite(psi).all() or not np.isfinite(g).all():
            raise ValueError('Finite corner psi and cell g required')
        scale = max(float(np.max(abs(psi))), np.finfo(float).tiny)
        if np.max(abs(psi[0])) > 1e-13*scale or np.ptp(psi[-1]) > 1e-13*scale or np.max(abs(psi[:, 0]-psi[:, -1])) > 1e-13*scale:
            raise ValueError('Regular axis, fixed wall and periodic flux required')
        old = geo.pack(psi, g).copy()
        ue, ui = model.electron.geometry.scalar(electron_entropy).copy(), model.ion.geometry.scalar(ion_entropy).copy()
        model.pressure(ue, ui)  # Validate before constructing any trial.
        jold = model.current(old)
        base_star = model.star_field(old, jold)
        impedance = max(float(np.max(abs(model.delta)*np.linalg.norm(base_star, axis=-1)/(model.density*E))),
                        2*float(model.inertia_coefficient.max())/dt)
        if impedance == 0:
            impedance = 1.  # Only a units/conditioning scale; no physical transport added.
        upper = (.5*dt*impedance*self.diagonal @ self.stiffness @ self.diagonal).tocsc()
        increment = np.zeros_like(old)
        candidate, ce, ci = old.copy(), ue.copy(), ui.copy()
        escale, iscale = float(abs(ue).max()), float(abs(ui).max())
        linear_error = 0.
        for iteration in range(1, self.max_iterations+1):
            mid = old+.5*increment
            jm = jold+model.current(.5*increment)
            em, im = logarithmic_mean(ue, ce), logarithmic_mean(ui, ci)
            pressure, _ = model.pressure(em, im)
            transport = model.delta*transport_matrix(geo, model.star_field(mid, jm), model.density)
            lower = -(self.diagonal @ (transport+2/dt*model.inertia) @ self.diagonal)/impedance
            matrix = bmat([[self.scaled_mass, upper], [lower, self.scaled_mass]], format='csc')
            rhs = np.r_[np.zeros_like(jold), self.scale*(transport @ jold+self.mass @ pressure)/impedance]
            solved = self._solve(matrix, rhs)
            linear_error = max(linear_error, float(np.linalg.norm(matrix @ solved-rhs)/max(np.linalg.norm(rhs), np.finfo(float).tiny)))
            electric = impedance*self.scale*solved[len(jold):]
            increment = dt*(geo.faraday @ electric)
            nextb = old+increment
            nextmid = old+.5*increment
            du, di = model.thermal_rates(geo.unpack(nextmid)[1], em, im)
            nexte, nexti = ue+dt*du, ui+dt*di
            if not np.isfinite(nextb).all() or not np.isfinite(nexte).all() or not np.isfinite(nexti).all() or np.any(nexte <= 0) or np.any(nexti <= 0):
                raise HallConvergenceError('Nonpositive/nonfinite current-thermal trial; no inputs committed')
            candidate, ce, ci = nextb, nexte, nexti
            mid = nextmid
            dj = model.current(.5*increment)
            jm = jold+dj
            em, im = logarithmic_mean(ue, ce), logarithmic_mean(ui, ci)
            pressure, _ = model.pressure(em, im)
            actual = model.delta*transport_matrix(geo, model.star_field(mid, jm), model.density)
            inertia = 2/dt*(model.inertia @ dj)
            terms = (self.mass @ electric, inertia, actual @ jm, self.mass @ pressure)
            ohm = terms[0]-terms[1]-terms[2]-terms[3]
            denominator = sum(np.linalg.norm(self.scale*v) for v in terms)
            ohm_error = float(np.linalg.norm(self.scale*ohm)/max(denominator, np.finfo(float).tiny))
            du, di = model.thermal_rates(geo.unpack(mid)[1], em, im)
            error = max(ohm_error, float(abs(ce-ue-dt*du).max())/escale, float(abs(ci-ui-dt*di).max())/iscale)
            if error <= self.tolerance and linear_error <= max(self.tolerance, 1e-12):
                # Check the stored state separately, on the matrix-product scale.
                # This exposes rounding when a small increment meets a large bias;
                # it is not a replacement for the small-signal Ohm residual above.
                stored_mid = .5*old+.5*candidate
                ampere = self.mass @ jm+geo.faraday.T @ (self.hessian @ stored_mid)
                envelope = abs(self.mass) @ abs(jm)+abs(geo.faraday.T) @ (abs(self.hessian) @ abs(stored_mid))
                ampere_error = np.linalg.norm(self.scale*ampere)/max(np.linalg.norm(self.scale*envelope), np.finfo(float).tiny)
                faraday_error = np.linalg.norm(candidate-old-increment)/max(np.linalg.norm(abs(candidate)+abs(old)+abs(increment)), np.finfo(float).tiny)
                stored_error = float(max(ampere_error, faraday_error))
                if not np.isfinite(stored_error) or stored_error > self.tolerance:
                    raise HallConvergenceError('Stored current-thermal state failed Ampere/Faraday backward residual; inputs unchanged')
                full_electric, powers = flux(model, mid, jm, em, im, electric)
                change = -dt*(np.diff(powers[0], axis=0)+np.diff(powers[1], axis=1))/mesh.volumes
                final_psi, final_g = geo.unpack(candidate)
                return CurrentThermalStep(final_psi, final_g, ce, ci, change, full_electric, *powers, iteration, error, linear_error, stored_error)
        raise HallConvergenceError(f'Joint current-thermal trial failed after {self.max_iterations} iterations; residual {error:.6g}; inputs unchanged')
