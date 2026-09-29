"""Unsplit conservative MHD in a fixed chamber coupled to a genuine vacuum."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from .coupling_mesh import CouplingMesh
from .mhd_state import MHDState, primitive
from .mhd_flux import fast_speed
from .vacuum_flux import flux_components
from .flux_limiting import admissible_update
from .vacuum_trace import tangential_field, interface_resistivity
from .vacuum_energy import variational_trace
from .magnetic_quadrature import selected_energy


@dataclass(frozen=True)
class ConservativeParameters:
    gamma: float = 5/3
    resistivity: float = 0.
    cfl: float = .2
    max_timestep: float | None = None
    vacuum_boundary: str = "multipole"
    multipole_modes: int = 12
    interface_diffusion: bool = True
    spatial_order: int = 2
    positivity_flux_limiting: bool = True
    flux_limiter: str = "local"
    energy_quadrature: str = "cell"
    boundary_poynting: str = "riemann"
    time_order: int = 2

    def __post_init__(self):
        if not np.isfinite(self.gamma) or self.gamma <= 1:
            raise ValueError("gamma must be finite and greater than one")
        if not np.isfinite(self.resistivity) or self.resistivity < 0:
            raise ValueError("Physical resistivity must be finite and nonnegative")
        if not np.isfinite(self.cfl) or not 0 < self.cfl <= .3:
            raise ValueError("Use 0 < CFL <= 0.3")
        if self.max_timestep is not None and (not np.isfinite(self.max_timestep) or self.max_timestep <= 0):
            raise ValueError("Source-resolution time cap must be finite and positive")
        if not isinstance(self.interface_diffusion, bool):
            raise ValueError("Interface numerical-diffusion switch must be boolean")
        if isinstance(self.spatial_order, bool) or self.spatial_order not in (1, 2):
            raise ValueError("Spatial order must be 1 or 2")
        if not isinstance(self.positivity_flux_limiting, bool):
            raise ValueError("Flux admissibility switch must be boolean")
        if self.flux_limiter not in ("local", "global"):
            raise ValueError("Flux limiter must be local or global")
        if self.energy_quadrature not in ("cell", "staggered"):
            raise ValueError("Energy quadrature must be cell or staggered")
        if self.boundary_poynting not in ("riemann", "variational"):
            raise ValueError("Boundary Poynting flux must be riemann or variational")
        if isinstance(self.time_order, bool) or self.time_order not in (2, 3):
            raise ValueError("Time order must be 2 or 3")


@dataclass
class ConservativeState:
    conserved: np.ndarray  # rho, rho vr, rho vz, total fluid-cell MHD energy
    induced_psi: np.ndarray  # full induced magnetic field grid, including vacuum

    def copy(self):
        return ConservativeState(self.conserved.copy(), self.induced_psi.copy())


class ConservativeVacuumMHD:
    """Reflecting gas chamber, transparent magnetic boundaries, no fluid exterior.

    Total-energy flux conservation is exact in the fluid cells. The independent
    full-space response-energy budget has measurable discretization error;
    it is never imposed through a thermal-energy repair.
    """
    def __init__(self, mesh, sources, initial, parameters=None):
        self.mesh, self.sources = mesh, sources
        self.parameters = parameters or ConservativeParameters()
        p = self.parameters
        if any(coil.radius <= mesh.radius and abs(coil.z) <= mesh.half_length for coil in sources.coils):
            raise ValueError("Prescribed source filaments must be outside the gas chamber")
        self.mapping = CouplingMesh(mesh, sources.grid, boundary=p.vacuum_boundary, multipole_modes=p.multipole_modes)
        self.state = initial.copy()
        self.time, self.steps = 0., 0
        self.euler_stages = 0
        self.source_work = self.ohmic_heat = self.outward_energy = 0.
        self.numerical_interface_heat = 0.
        self.limited_euler_stages, self.minimum_high_order_weight = 0, 1.
        self.limiter_cell_counts = np.zeros((mesh.nr, mesh.nz), dtype=np.int64)
        self.limiter_troubled_counts = np.zeros_like(self.limiter_cell_counts)
        self.maximum_limiter_cell_fraction = 0.
        self.limiter_iterations = self.maximum_limiter_iterations = 0
        self._primitive(self.state, 0.)
        u = self.mapping.interior(initial.induced_psi)
        ku = self.mapping.k@u
        if np.linalg.norm(ku[self.mapping.vacuum]) > 1e-10*max(np.linalg.norm(ku), 1e-30):
            raise ValueError("Initial induced field violates the vacuum constraint")
        if np.any(initial.induced_psi[0]):
            raise ValueError("Axis induced flux must vanish")
        if not self.mapping.open_boundary and (np.any(initial.induced_psi[-1]) or np.any(initial.induced_psi[:, [0, -1]])):
            raise ValueError("Fixed outer induced flux must vanish")
        self.initial = self._integrals(self.state, 0.)

    def total_flux(self, state=None, time=None):
        state = self.state if state is None else state
        time = self.time if time is None else time
        return state.induced_psi+self.sources.at(time).psi

    def _primitive(self, state, time):
        flux = self.mapping.fluid_flux(self.total_flux(state, time))
        return primitive(self.mesh, MHDState(state.conserved, flux), self.parameters.gamma,
                         energy_quadrature=self.parameters.energy_quadrature)

    def _integrals(self, state, time):
        q = self._primitive(state, time)
        vol, m = self.mesh.volumes, self.mapping
        flux = m.fluid_flux(self.total_flux(state, time))
        magnetic_cell = float(np.sum(selected_energy(self.mesh, flux, self.parameters.energy_quadrature)*vol))
        energy = float(np.sum(state.conserved[..., 3]*vol))
        u = m.interior(state.induced_psi)
        induced_energy = float(np.pi/MU0*(u@(m.k@u)))
        return {"mass_kg": float(np.sum(q[..., 0]*vol)), "fluid_total_energy_j": energy,
                "gas_energy_j": energy-magnetic_cell, "fluid_magnetic_energy_j": magnetic_cell,
                "induced_field_energy_j": induced_energy, "response_energy_j": energy-magnetic_cell+induced_energy,
                "internal_energy_j": float(np.sum(q[..., 3]*vol/(self.parameters.gamma-1))),
                "minimum_density_kg_m3": float(q[..., 0].min()), "minimum_pressure_pa": float(q[..., 3].min()),
                "maximum_speed_m_s": float(np.hypot(q[..., 1], q[..., 2]).max())}

    def timestep_limit(self):
        q = self._primitive(self.state, self.time)
        mesh, p = self.mesh, self.parameters
        waves = np.max((abs(q[..., 1])+fast_speed(q, 0, p.gamma))/mesh.dr+
                       (abs(q[..., 2])+fast_speed(q, 1, p.gamma))/mesh.dz)
        eta = interface_resistivity(self.mapping, q, p.gamma) if p.interface_diffusion else np.zeros(1)
        diffusion = 2*(p.resistivity+float(eta.max()))/MU0*(1/mesh.dr**2+1/mesh.dz**2)
        limit = float(p.cfl/(waves+diffusion))
        return limit if p.max_timestep is None else min(limit, p.max_timestep)

    def _fluxes(self, state, time, spatial_order=None):
        m, p = self.mapping, self.parameters
        q = self._primitive(state, time)
        ku = m.k@m.interior(state.induced_psi)
        current = np.divide(ku, MU0*m.node_area, out=np.zeros(m.node_count), where=m.active)
        electric = m.fluid_flux(m.full(p.resistivity*current))
        full_total = self.total_flux(state, time)
        total = m.fluid_flux(full_total)
        numerical_eta = interface_resistivity(m, q, p.gamma) if p.interface_diffusion else np.zeros(m.node_count)
        numerical_electric = m.fluid_flux(m.full(numerical_eta*current))
        boundary = variational_trace(m, full_total, ku) if p.boundary_poynting == "variational" else None
        traces = boundary.face_fields(self.mesh) if boundary is not None else tangential_field(m, full_total)
        fluxes = flux_components(self.mesh, q, total, p.gamma, electric, traces,
                                  numerical_electric, p.spatial_order if spatial_order is None else spatial_order)
        fluxes.boundary_poynting = boundary
        source_derivative = m.grid_vector(self.sources.at(time).dpsi_dt)
        source_power = -float((2*np.pi/MU0)*ku[m.active]@source_derivative[m.active])
        ohmic_power = float(2*np.pi*np.sum(p.resistivity*current**2*m.node_radius*m.node_area))
        numerical_power = float(2*np.pi*np.sum(numerical_eta*current**2*m.node_radius*m.node_area))
        return fluxes, (source_power, ohmic_power, numerical_power)

    def _forward_euler(self, state, time, total, dt):
        """Select shared face/corner fluxes before assembly; never repair energy."""
        high, powers = self._fluxes(state, time)
        p = self.parameters
        mode = p.flux_limiter if p.positivity_flux_limiting and p.spatial_order == 2 else "off"
        rhs, report = admissible_update(self.mesh, state.conserved, total, dt, high,
                                        lambda: self._fluxes(state, time, spatial_order=1)[0], mode=mode,
                                        energy_quadrature=p.energy_quadrature)
        # Current, physical/interface eta and analytic source power are independent
        # of reconstruction order at this state; only outward flux is selected.
        return (*rhs, *powers), report

    def _new_state(self, conserved, fluid_total_flux, time):
        m, mesh = self.mapping, self.mesh
        data = np.zeros(m.field.shape)
        external = self.sources.at(time).psi
        data[:mesh.nr+1, m.z_offset:m.z_offset+mesh.nz+1] = fluid_total_flux-m.fluid_flux(external)
        result = ConservativeState(conserved, m.harmonic_extension(data))
        self._primitive(result, time)
        return result

    def step(self, dt=None):
        limit = self.timestep_limit()
        dt = limit if dt is None else dt
        if not np.isfinite(dt) or dt <= 0 or dt > limit*(1+1e-12):
            raise ValueError(f"Time step must be positive and <= {limit:.6g} s")
        old = self.state
        total = self.mapping.fluid_flux(self.total_flux())
        a, report_a = self._forward_euler(old, self.time, total, dt)
        predictor = self._new_state(old.conserved+dt*a[0], total+dt*a[1], self.time+dt)
        b, report_b = self._forward_euler(predictor, self.time+dt, total+dt*a[1], dt)
        if self.parameters.time_order == 2:
            final = self._new_state(old.conserved+.5*dt*(a[0]+b[0]), total+.5*dt*(a[1]+b[1]), self.time+dt)
            weighted = ((.5, a), (.5, b))
            reports = (report_a, report_b)
        else:
            half_flux = total+.25*dt*(a[1]+b[1])
            half = self._new_state(old.conserved+.25*dt*(a[0]+b[0]), half_flux, self.time+.5*dt)
            c, report_c = self._forward_euler(half, self.time+.5*dt, half_flux, dt)
            final = self._new_state(old.conserved+dt*((a[0]+b[0])/6+2*c[0]/3),
                                    total+dt*((a[1]+b[1])/6+2*c[1]/3), self.time+dt)
            weighted = ((1/6, a), (1/6, b), (2/3, c))
            reports = (report_a, report_b, report_c)
        # Every Euler proposal and exact vacuum completion validates atomically.
        self.state = final
        self.time += dt; self.steps += 1
        self.euler_stages += len(reports)
        self.outward_energy += dt*sum(w*rates[2] for w, rates in weighted)
        self.source_work += dt*sum(w*rates[3] for w, rates in weighted)
        self.ohmic_heat += dt*sum(w*rates[4] for w, rates in weighted)
        self.numerical_interface_heat += dt*sum(w*rates[5] for w, rates in weighted)
        for report in reports:
            self.limited_euler_stages += int(report.minimum_weight < 1)
            self.minimum_high_order_weight = min(self.minimum_high_order_weight, report.minimum_weight)
            self.limiter_cell_counts += report.affected
            self.limiter_troubled_counts += report.troubled
            self.maximum_limiter_cell_fraction = max(self.maximum_limiter_cell_fraction, float(report.affected.mean()))
            self.limiter_iterations += report.iterations
            self.maximum_limiter_iterations = max(self.maximum_limiter_iterations, report.iterations)
        return dt

    def advance(self, end_time, *, max_steps=1_000_000, on_step=None):
        if not np.isfinite(end_time) or end_time < self.time:
            raise ValueError("End time must be finite and not precede accepted time")
        count = 0
        while self.time < end_time:
            remaining = end_time-self.time
            if remaining <= 4*np.finfo(float).eps*max(abs(end_time), 1e-300):
                self.time = end_time; break
            if count >= max_steps:
                raise RuntimeError("Step budget exceeded; accepted state retained")
            self.step(min(self.timestep_limit(), remaining)); count += 1
            if on_step is not None:
                on_step(self)
        return self.state

    def diagnostics(self):
        result = self._integrals(self.state, self.time)
        response_residual = result["response_energy_j"]-self.initial["response_energy_j"]-self.source_work
        fluid_residual = result["fluid_total_energy_j"]-self.initial["fluid_total_energy_j"]+self.outward_energy
        ku = self.mapping.k@self.mapping.interior(self.state.induced_psi)
        result.update(time_s=self.time, steps=self.steps, euler_stages=self.euler_stages,
                      source_work_j=self.source_work, ohmic_heat_j=self.ohmic_heat,
                      numerical_interface_heat_j=self.numerical_interface_heat,
                      limited_euler_stages=self.limited_euler_stages, minimum_high_order_weight=self.minimum_high_order_weight,
                      affected_euler_cell_stages=int(self.limiter_cell_counts.sum()),
                      troubled_euler_cell_stages=int(self.limiter_troubled_counts.sum()),
                      maximum_limiter_cell_fraction=self.maximum_limiter_cell_fraction,
                      limiter_iterations=self.limiter_iterations, maximum_limiter_iterations=self.maximum_limiter_iterations,
                      outward_fluid_energy_j=self.outward_energy, fluid_energy_balance_residual_j=fluid_residual,
                      response_energy_balance_residual_j=response_residual,
                      response_residual_over_source_work=response_residual/abs(self.source_work) if self.source_work else None,
                      response_energy_balance_relative=response_residual/max(abs(self.initial["response_energy_j"]), abs(self.source_work), 1e-30),
                      mass_relative_change=(result["mass_kg"]-self.initial["mass_kg"])/self.initial["mass_kg"],
                      vacuum_relative_residual=float(np.linalg.norm(ku[self.mapping.vacuum])/max(np.linalg.norm(ku), 1e-30)),
                      max_face_divergence_t_m=float(abs(self.mesh.divergence(self.mapping.fluid_flux(self.total_flux()))).max()))
        return result
