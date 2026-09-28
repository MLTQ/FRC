"""Developmental gas dynamics coupled to prescribed coils and true vacuum."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from .coupling_mesh import CouplingMesh
from .coupled_hydro import gas_primitive, hydro_step
from .coupled_magnetic import MagneticStepper
from .numerics import cell_divergence
from .magnetic_dissipation import shock_resistance


@dataclass(frozen=True)
class CoupledParameters:
    gamma: float = 5/3
    resistivity: float = 0.
    cfl: float = .2
    nonlinear_tolerance: float = 1e-13
    max_nonlinear_iterations: int = 40
    max_timestep: float | None = None
    vacuum_boundary: str = "fixed"
    multipole_modes: int = 12
    midpoint_method: str = "picard"
    shock_diffusion: float = 0.

    def __post_init__(self):
        if not np.isfinite(self.gamma) or self.gamma <= 1:
            raise ValueError("gamma must be finite and greater than one")
        if not np.isfinite(self.cfl) or not 0 < self.cfl <= .3:
            raise ValueError("Use 0 < CFL <= 0.3 for the coupled reference")
        if self.max_timestep is not None and (not np.isfinite(self.max_timestep) or self.max_timestep <= 0):
            raise ValueError("An optional source-resolution time cap must be finite and positive")
        if self.vacuum_boundary not in ("fixed", "multipole"):
            raise ValueError("Magnetic boundary must be 'fixed' or 'multipole'")
        if self.midpoint_method not in ("direct", "picard"):
            raise ValueError("Midpoint method must be 'direct' or 'picard'")
        if not np.isfinite(self.shock_diffusion) or not 0 <= self.shock_diffusion <= 2:
            raise ValueError("Use a declared numerical shock-diffusion coefficient between zero and two")


@dataclass
class CoupledState:
    gas: np.ndarray  # rho, rho vr, rho vz, GAS total energy (no magnetic energy)
    induced_psi: np.ndarray  # full plasma+vacuum field grid

    def copy(self):
        return CoupledState(self.gas.copy(), self.induced_psi.copy())


class VacuumCoupledMHD:
    """Fixed gas-filled chamber, reflecting ends, magnetically transparent walls.

    An evolving detached-plasma interface and physical transport are NOT yet
    implemented; the original Stage 2 formation entry point remains gated.
    """
    def __init__(self, fluid_mesh, vacuum_field, initial: CoupledState, parameters=None):
        self.mesh, self.sources = fluid_mesh, vacuum_field
        self.parameters = parameters or CoupledParameters()
        if any(coil.radius <= fluid_mesh.radius and abs(coil.z) <= fluid_mesh.half_length
               for coil in vacuum_field.coils):
            raise ValueError("Prescribed source filaments must be outside the gas-filled chamber")
        self.mapping = CouplingMesh(fluid_mesh, vacuum_field.grid, boundary=self.parameters.vacuum_boundary,
                                    multipole_modes=self.parameters.multipole_modes)
        self.magnetic = MagneticStepper(self.mapping, self.parameters.resistivity, self.parameters.gamma,
                                        tolerance=self.parameters.nonlinear_tolerance,
                                        max_iterations=self.parameters.max_nonlinear_iterations,
                                        method=self.parameters.midpoint_method, shock_diffusion=self.parameters.shock_diffusion)
        gas_primitive(self.mesh, initial.gas, self.parameters.gamma)
        u = self.mapping.interior(initial.induced_psi)
        if np.any(initial.induced_psi[0]):
            raise ValueError("Induced flux must vanish on the symmetry axis")
        if not self.mapping.open_boundary and (np.any(initial.induced_psi[-1]) or np.any(initial.induced_psi[:, [0, -1]])):
            raise ValueError("Induced flux must vanish on fixed outer field boundaries")
        ku = self.mapping.k@u
        if np.linalg.norm(ku[self.mapping.vacuum]) > 1e-10*max(np.linalg.norm(ku), 1e-30):
            raise ValueError("Initial induced flux violates the exterior vacuum constraint")
        self.state = initial.copy()
        self.time, self.steps = 0., 0
        self.source_work = self.ohmic_heat = self.lorentz_work = self.exchange_residual = 0.
        self.numerical_heat = self.maximum_numerical_resistivity = 0.
        self.max_iterations = 0
        self.max_nonlinear_residual = 0.
        self.initial = self._integrals(self.state)

    def total_flux(self, state=None, time=None):
        state = self.state if state is None else state
        return state.induced_psi+self.sources.at(self.time if time is None else time).psi

    def _integrals(self, state):
        q = gas_primitive(self.mesh, state.gas, self.parameters.gamma)
        volume = self.mesh.volumes
        u = self.mapping.interior(state.induced_psi)
        gas_energy = float(np.sum(state.gas[..., 3]*volume))
        field_energy = float(np.pi/MU0*(u@(self.mapping.k@u)))
        exterior_energy = 0.
        if self.mapping.open_boundary:
            modes = u[self.mapping.grid_node_count:]
            exterior_energy = float(np.pi/MU0*(modes@self.mapping.open_boundary.exterior@modes))
        return {"mass_kg": float(np.sum(q[..., 0]*volume)), "gas_energy_j": gas_energy,
                "internal_energy_j": float(np.sum(q[..., 3]*volume/(self.parameters.gamma-1))),
                "kinetic_energy_j": float(np.sum(.5*q[..., 0]*(q[..., 1]**2+q[..., 2]**2)*volume)),
                "induced_field_energy_j": field_energy, "response_energy_j": gas_energy+field_energy,
                "beyond_grid_induced_energy_j": exterior_energy,
                "minimum_density_kg_m3": float(q[..., 0].min()), "minimum_pressure_pa": float(q[..., 3].min()),
                "maximum_speed_m_s": float(np.hypot(q[..., 1], q[..., 2]).max())}

    def timestep_limit(self):
        q = gas_primitive(self.mesh, self.state.gas, self.parameters.gamma)
        br, bz = self.mesh.cell_fields(self.mapping.fluid_flux(self.total_flux()))
        speed = np.sqrt(self.parameters.gamma*q[..., 3]/q[..., 0]+(br*br+bz*bz)/(MU0*q[..., 0]))
        waves = np.max((abs(q[..., 1])+speed)/self.mesh.dr+(abs(q[..., 2])+speed)/self.mesh.dz)
        # Midpoint diffusion is implicit, but this accuracy bound avoids poorly
        # resolved alternating high-frequency diffusion modes.
        maximum_eta = 0.
        if self.parameters.shock_diffusion:
            velocity = np.concatenate((q[..., 1].ravel(), q[..., 2].ravel()))
            _, maximum_eta = shock_resistance(self.mapping, q, self.total_flux(), velocity,
                                               self.parameters.shock_diffusion, self.parameters.gamma)
        diffusion = 2*(self.parameters.resistivity+maximum_eta)/MU0*(1/self.mesh.dr**2+1/self.mesh.dz**2)
        limit = float(self.parameters.cfl/(waves+diffusion))
        return limit if self.parameters.max_timestep is None else min(limit, self.parameters.max_timestep)

    def step(self, dt=None):
        limit = self.timestep_limit()
        dt = limit if dt is None else dt
        if not np.isfinite(dt) or dt <= 0 or dt > limit*(1+1e-12):
            raise ValueError(f"Coupled time step must be positive and <= {limit:.6g} s")
        first = hydro_step(self.mesh, self.state.gas, dt/2, self.parameters.gamma)
        update = self.magnetic.step(first, self.state.induced_psi, self.sources.at(self.time).psi,
                                    self.sources.at(self.time+dt).psi, dt)
        final = hydro_step(self.mesh, update.gas, dt/2, self.parameters.gamma)
        # Commit only after all three split stages have validated.
        self.state = CoupledState(final, update.induced)
        self.time += dt; self.steps += 1
        self.source_work += update.source_work
        self.ohmic_heat += update.ohmic_heat
        self.lorentz_work += update.lorentz_work
        self.numerical_heat += update.numerical_heat
        self.maximum_numerical_resistivity = max(self.maximum_numerical_resistivity, update.maximum_numerical_resistivity)
        self.exchange_residual += update.exchange_residual
        self.max_iterations = max(self.max_iterations, update.iterations)
        self.max_nonlinear_residual = max(self.max_nonlinear_residual, update.nonlinear_residual)
        return dt

    def advance(self, end_time, *, max_steps=1_000_000):
        if not np.isfinite(end_time) or end_time < self.time:
            raise ValueError("End time must be finite and not precede accepted state")
        count = 0
        while self.time < end_time:
            remaining = end_time-self.time
            if remaining <= 4*np.finfo(float).eps*max(abs(end_time), 1e-300):
                self.time = end_time; break
            if count >= max_steps:
                raise RuntimeError("Coupled step budget exceeded; accepted state is retained")
            self.step(min(self.timestep_limit(), remaining)); count += 1
        return self.state

    def diagnostics(self):
        out = self._integrals(self.state)
        ku = self.mapping.k@self.mapping.interior(self.state.induced_psi)
        grid_nodes = self.mapping.grid_node_count
        vacuum_j = ku[:grid_nodes][self.mapping.vacuum[:grid_nodes]]/(MU0*self.mapping.field.dr*self.mapping.field.dz)
        exterior_residual = float(np.linalg.norm(ku[grid_nodes:])/max(np.linalg.norm(ku[self.mapping.active]), 1e-30))
        residual = out["response_energy_j"]-self.initial["response_energy_j"]-self.source_work
        out.update(time_s=self.time, steps=self.steps, source_work_j=self.source_work,
                   ohmic_heat_j=self.ohmic_heat, lorentz_work_j=self.lorentz_work,
                   numerical_magnetic_heat_j=self.numerical_heat,
                   maximum_numerical_resistivity_ohm_m=self.maximum_numerical_resistivity,
                   response_energy_balance_residual_j=residual,
                   response_energy_balance_relative=residual/max(abs(self.initial["response_energy_j"]), abs(self.source_work), 1e-30),
                   mass_relative_change=(out["mass_kg"]-self.initial["mass_kg"])/self.initial["mass_kg"],
                   magnetic_exchange_residual_j=self.exchange_residual,
                   max_midpoint_iterations=self.max_iterations, max_midpoint_relative_residual=self.max_nonlinear_residual,
                   open_boundary_relative_stationarity_residual=exterior_residual,
                   max_vacuum_induced_current_residual_a_m2=float(abs(vacuum_j).max(initial=0)),
                   max_face_divergence_t_m=float(abs(cell_divergence(self.mapping.field, self.total_flux())).max()))
        return out
