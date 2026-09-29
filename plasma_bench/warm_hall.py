"""Warm massless-electron Hall MHD with explicit conservative thermal transport."""
import numpy as np
from scipy.constants import elementary_charge as E
from .warm_mhd import AdiabaticSpeciesMHD
from .hall_trace import TraceHallCurl
from .hall_midpoint import HallMidpoint
from .electron_pressure import ElectronPressure
from .electron_midpoint import step as pressure_step
from .electron_energy import energy_rate as pressure_energy_rate
from .thermal_means import logarithmic_mean
from .compatible_poynting import energy_rate as hall_energy_rate


class WarmHallMHD(AdiabaticSpeciesMHD):
    """Ideal, equal-gamma, positive filled cylinder; fixed psi and zero wall Btheta."""
    def __init__(self, mesh, initial, electron_pressure, ion_mass, parameters=None, *, density_scale=1e20,
                 quadrature_order=5, nonlinear_tolerance=1e-10, hall_iterations=16, pressure_iterations=80):
        if initial.magnetic_quadrature != 'compatible':
            raise ValueError('Warm Hall requires explicit compatible magnetic energy')
        if type(pressure_iterations) is not int or pressure_iterations < 1:
            raise ValueError('Positive pressure iteration budget required')
        super().__init__(mesh, initial, electron_pressure, ion_mass, parameters, density_scale=density_scale)
        self.geometry = TraceHallCurl(mesh, quadrature_order)
        self.induction = HallMidpoint(self.geometry, tolerance=nonlinear_tolerance, max_iterations=hall_iterations)
        self.pressure_tolerance, self.pressure_budget = nonlinear_tolerance, pressure_iterations
        self.hall_iterations = self.pressure_iterations = 0
        self.max_hall_residual = self.max_linear_residual = self.max_pressure_residual = 0.

    def species_rhs(self, state, entropy):
        derivative, du, heat, statistics = super().species_rhs(state, entropy)
        density = self.thermal.density(state)
        model = ElectronPressure(self.mesh, density, gamma=self.parameters.gamma, density_scale=self.thermal.density_scale)
        dg, drift, _, _ = model.rhs(state.conserved[..., 5], entropy)
        current = self.geometry.current(state.psi, state.conserved[..., 5])
        magnetic = self.geometry.magnetic_samples(state.psi, state.conserved[..., 5])
        electric = np.cross(current, magnetic)/(E*self.geometry.scalar_samples(density)[..., None])
        coefficients = self.geometry.project_electric(electric)
        dp, hall_g = self.geometry.unpack(self.geometry.faraday @ coefficients)
        derivative.psi += dp
        derivative.conserved[..., 5] += hall_g+dg
        derivative.conserved[..., 3] += hall_energy_rate(self.geometry, state.psi, state.conserved[..., 5], coefficients)
        derivative.conserved[..., 3] += pressure_energy_rate(model, state.conserved[..., 5], entropy)
        return derivative, du+drift, heat, statistics

    def rhs(self, state, entropy_density=None):
        if entropy_density is None:
            raise ValueError('Warm Hall RHS requires the matching electron entropy state')
        derivative, _, heat, _ = self.species_rhs(state, entropy_density)
        return derivative, heat

    def _ideal(self, state, entropy, dt):
        trial = AdiabaticSpeciesMHD(self.mesh, state, self.thermal.pressure(state, entropy), self.thermal.ion_mass,
                                    self.parameters, density_scale=self.thermal.density_scale)
        trial.electron_entropy = entropy.copy()
        trial.step(dt)
        return trial.state, trial.electron_entropy, trial.limiter_totals

    def _pressure(self, state, entropy, dt):
        model = ElectronPressure(self.mesh, self.thermal.density(state), gamma=self.parameters.gamma,
                                  density_scale=self.thermal.density_scale, wall_btheta=0.)
        oldg = state.conserved[..., 5]
        result = pressure_step(model, oldg, entropy, dt, tolerance=self.pressure_tolerance, max_iterations=self.pressure_budget)
        gmid = .5*(oldg+result.g)
        umid = logarithmic_mean(entropy, result.entropy_density)
        trial = state.copy()
        trial.conserved[..., 5] = result.g
        trial.conserved[..., 3] += dt*pressure_energy_rate(model, gmid, umid)
        self.thermal.validate(trial, result.entropy_density)
        return trial, result

    def step(self, dt=None):
        limit = self.timestep_limit()
        dt = limit if dt is None else dt
        if not np.isfinite(dt) or dt <= 0 or dt > limit*(1+1e-12):
            raise ValueError(f'Positive time step <= fluid limit {limit:.6g} required')
        first, entropy, stats_a = self._ideal(self.state, self.electron_entropy, .5*dt)
        first, pressure_a = self._pressure(first, entropy, .5*dt)
        entropy = pressure_a.entropy_density
        density = self.geometry.scalar_samples(self.thermal.density(first))
        result = self.induction.step(first.psi, first.conserved[..., 5], density, dt)
        middle = first.copy()
        middle.psi = result.psi
        middle.conserved[..., 5] = result.g
        middle.conserved[..., 3] += result.total_energy_change
        self.thermal.validate(middle, entropy)
        middle, pressure_b = self._pressure(middle, entropy, .5*dt)
        final, entropy, stats_b = self._ideal(middle, pressure_b.entropy_density, .5*dt)
        # No accepted state/account changes until all five substeps succeed.
        self.state, self.electron_entropy = final, entropy
        self.time += dt; self.steps += 1
        self.hall_iterations += result.iterations
        self.pressure_iterations += pressure_a.iterations+pressure_b.iterations
        self.max_hall_residual = max(self.max_hall_residual, result.nonlinear_residual)
        self.max_linear_residual = max(self.max_linear_residual, result.linear_residual)
        self.max_pressure_residual = max(self.max_pressure_residual, pressure_a.residual, pressure_b.residual)
        for key in self.limiter_totals:
            if key == 'minimum_pressure_blend_weight':
                self.limiter_totals[key] = min(self.limiter_totals[key], stats_a[key], stats_b[key])
            else:
                self.limiter_totals[key] += stats_a[key]+stats_b[key]
        return dt

    def diagnostics(self):
        result = super().diagnostics()
        result.update(species_model='equal-gamma massless-electron Hall with pressure battery/enthalpy; eta=0',
                      magnetic_boundary='fixed psi and Btheta(R)=0; free Ez(R)',
                      hall_nonlinear_iterations=self.hall_iterations, pressure_nonlinear_iterations=self.pressure_iterations,
                      max_hall_nonlinear_residual=self.max_hall_residual, max_hall_linear_residual=self.max_linear_residual,
                      max_pressure_nonlinear_residual=self.max_pressure_residual)
        result['reconstruction_limiter']['accepted_cell_stages'] = 4*self.steps*self.mesh.nr*self.mesh.nz
        return result
