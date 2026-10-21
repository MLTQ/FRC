"""Unsplit conservative MHD in a fixed chamber coupled to a genuine vacuum."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0 as MU0
from .coupling_mesh import CouplingMesh
from .mhd_state import MHDState, primitive, toroidal_energies
from .mhd_flux import fast_speed
from .toroidal_flux import fast_speed as toroidal_fast_speed
from . import toroidal_ohm
from .vacuum_flux import flux_components
from .flux_limiting import admissible_update
from .vacuum_trace import tangential_field, interface_resistivity
from .vacuum_energy import variational_trace
from .magnetic_quadrature import selected_energy
from .dual_energy import entropy_density, entropy_rate, cold_cells, select
from .mhd_state import NonphysicalState
from .density_floor import apply_floor
from .open_line_loss import open_line_loss


@dataclass(frozen=True)
class ConservativeParameters:
    gamma: float = 5/3
    resistivity: float = 0.
    resistivity_model: str = "uniform"  # "uniform" (scalar above) or "spitzer" (NRL eta_perp of local Te)
    spitzer_coulomb_log: float = 10.
    spitzer_min_te_ev: float = 1.       # Te floor for evaluating Spitzer eta only (no state change)
    spitzer_ion_mass: float = 2.01410*1.66053906660e-27  # kg per ion (Z=1, Te=Ti)
    chodura_coefficient: float = 0.     # >0 adds Chodura anomalous eta = c m_e w_pi/(n e^2) (1-exp(-v_d/(f c_s)))
    chodura_f: float = 3.
    gap_resistivity: float = 0.         # >0: DECLARED resistivity [ohm m] of near-floor cells (un-ionized gas / gap):
    gap_density: float | None = None    # weight 1/(1+(rho/gap_density)^4) at nodes; makes them magnetically transparent
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
    dual_energy: bool = False           # opt-in entropy pressure in cold/low-beta cells
    dual_energy_threshold: float = .1   # internal/total energy below which a cell is cold
    dual_energy_mode: str = "pressure"  # entropy pressure only (energy exact) or "reset" total energy
    dual_energy_shock_threshold: float | None = .2  # compression detector keeps shocks on total energy
    dual_energy_shock_min_fraction: float = .01  # detector needs internal >= this * magnetic energy (low beta stays cold)
    dual_energy_criterion: str = "total"  # cold test on E-K-M ("total") or on the entropy thermal energy ("entropy")
    max_step_halvings: int = 0          # retry an inadmissible step from the accepted state with dt/2
    local_trace_bounds: bool = False    # reconstructed density/pressure >= 3x3 neighbour minimum
    interface_speed_cap: float | None = None  # bound on the signal speed in the interface numerical resistivity
    open_ends: bool = False             # outflow (zero-gradient, no inflow) gas ends instead of reflecting ends
    analytic_source_traces: bool = False  # wall/end traces: differenced induced flux + analytic coil field
    magnetic_energy_flux: str = "face"  # "dissipative": LLF magnetic dissipation as the Poynting flux of the
                                        # UCT dissipative field; "ct": all magnetic energy transport as the
                                        # Poynting flux of the full CT electric field (vacuum_flux.md)
    open_line_flux_limit: float = 0.    # >0: DECLARED flux-limited parallel conduction loss on open lines (f, ~0.1)
    open_line_wall_temperature_ev: float = 2.  # end-wall temperature the open-line plasma relaxes toward
    minimum_density: float | None = None  # DECLARED gap closure [kg/m^3]; added mass/energy are counted
    floor_temperature_ev: float = 1.     # temperature of particles added by the closure
    floor_ion_mass: float = 2.01410*1.66053906660e-27  # kg per ion for the closure's thermal energy
    toroidal: bool = False              # 2.5D: also evolve r rho v_theta and B_theta/r; insulating gas boundaries
                                        # (B_theta = 0 there unless `toroidal_boundary` threads a current)
    hall: bool = False                  # DECLARED Hall term E += J x B/(n e) (needs toroidal and time_order=3)
    hall_ion_mass: float | None = None  # kg per ion for n = rho/m_i in the Hall term (None: spitzer_ion_mass)
    hall_cutoff_density: float | None = None  # DECLARED Hall cutoff [kg/m^3]: weight 1/(1+(rho_c/rho)^4)
    hall_whistler_fraction: float = .5  # whistler step: this fraction of sqrt(3)/omega_max (RK3 imaginary axis)
    toroidal_boundary: object | None = None  # callable t -> r B_theta [T m] of the vacuum beyond every gas boundary
                                        # (mu0 I/(2 pi) for a current I threading an annular chamber); None: zero

    def __post_init__(self):
        if not np.isfinite(self.gamma) or self.gamma <= 1:
            raise ValueError("gamma must be finite and greater than one")
        if not np.isfinite(self.resistivity) or self.resistivity < 0:
            raise ValueError("Physical resistivity must be finite and nonnegative")
        if self.resistivity_model not in ("uniform", "spitzer"):
            raise ValueError("Resistivity model must be uniform or spitzer")
        if not (self.spitzer_coulomb_log > 0 and self.spitzer_min_te_ev > 0 and self.spitzer_ion_mass > 0):
            raise ValueError("Spitzer parameters must be positive")
        if not (np.isfinite(self.chodura_coefficient) and self.chodura_coefficient >= 0 and self.chodura_f > 0):
            raise ValueError("Chodura coefficient must be nonnegative and f positive")
        if not (np.isfinite(self.gap_resistivity) and self.gap_resistivity >= 0):
            raise ValueError("Gap resistivity must be finite and nonnegative")
        if self.gap_resistivity > 0 and not (self.gap_density is not None and self.gap_density > 0):
            raise ValueError("Gap resistivity needs a positive gap density")
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
        if self.magnetic_energy_flux not in ("face", "dissipative", "ct"):
            raise ValueError("Magnetic energy flux must be face, dissipative or ct")
        if self.boundary_poynting not in ("riemann", "variational"):
            raise ValueError("Boundary Poynting flux must be riemann or variational")
        if isinstance(self.time_order, bool) or self.time_order not in (2, 3):
            raise ValueError("Time order must be 2 or 3")
        if self.minimum_density is not None and not (np.isfinite(self.minimum_density) and self.minimum_density > 0):
            raise ValueError("Minimum density must be positive or None")
        if not (np.isfinite(self.floor_temperature_ev) and self.floor_temperature_ev > 0 and self.floor_ion_mass > 0):
            raise ValueError("Floor temperature and ion mass must be positive")
        if not (np.isfinite(self.open_line_flux_limit) and self.open_line_flux_limit >= 0 and self.open_line_wall_temperature_ev > 0):
            raise ValueError("Open-line flux limit must be nonnegative and wall temperature positive")
        if not isinstance(self.local_trace_bounds, bool):
            raise ValueError("Local trace-bound switch must be boolean")
        if not isinstance(self.dual_energy, bool):
            raise ValueError("Dual-energy switch must be boolean")
        if not np.isfinite(self.dual_energy_threshold) or not 0 < self.dual_energy_threshold < .5:
            raise ValueError("Use 0 < dual-energy threshold < 0.5")
        if self.dual_energy_criterion not in ("total", "entropy"):
            raise ValueError("Dual-energy criterion must be total or entropy")
        if self.dual_energy_mode not in ("reset", "pressure"):
            raise ValueError("Dual-energy mode must be reset or pressure")
        if isinstance(self.max_step_halvings, bool) or not isinstance(self.max_step_halvings, int) or self.max_step_halvings < 0:
            raise ValueError("Step halvings must be a nonnegative integer")
        if self.dual_energy_shock_threshold is not None and not (np.isfinite(self.dual_energy_shock_threshold)
                                                                 and self.dual_energy_shock_threshold > 0):
            raise ValueError("Dual-energy shock threshold must be positive or None")
        if not isinstance(self.toroidal, bool) or not isinstance(self.hall, bool):
            raise ValueError("Toroidal and Hall switches must be boolean")
        if self.hall and not self.toroidal:
            raise ValueError("The Hall term generates B_theta and J_r, J_z: it needs toroidal=True")
        if self.hall and self.time_order != 3:
            raise ValueError("The Hall term needs time_order=3: RK2 amplifies undamped whistlers")
        if self.hall_ion_mass is not None and not (np.isfinite(self.hall_ion_mass) and self.hall_ion_mass > 0):
            raise ValueError("Hall ion mass must be positive or None")
        if self.hall_cutoff_density is not None and not (np.isfinite(self.hall_cutoff_density) and self.hall_cutoff_density > 0):
            raise ValueError("Hall cutoff density must be positive or None")
        if not (np.isfinite(self.hall_whistler_fraction) and 0 < self.hall_whistler_fraction <= 1):
            raise ValueError("Use 0 < Hall whistler fraction <= 1")
        if self.toroidal_boundary is not None and not (callable(self.toroidal_boundary) and self.toroidal):
            raise ValueError("A boundary r B_theta must be a callable of time and needs toroidal=True")


@dataclass
class ConservativeState:
    conserved: np.ndarray  # rho, rho vr, rho vz, total fluid-cell MHD energy [, L = r rho v_theta, g = B_theta/r]
    induced_psi: np.ndarray  # full induced magnetic field grid, including vacuum
    entropy: np.ndarray | None = None  # rho*p/rho^gamma; only with dual energy

    def copy(self):
        return ConservativeState(self.conserved.copy(), self.induced_psi.copy(),
                                 None if self.entropy is None else self.entropy.copy())


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
        if any(mesh.inner_radius <= coil.radius <= mesh.radius and abs(coil.z) <= mesh.half_length for coil in sources.coils):
            raise ValueError("Prescribed source filaments must be outside the gas chamber")
        self.mapping = CouplingMesh(mesh, sources.grid, boundary=p.vacuum_boundary, multipole_modes=p.multipole_modes)
        self.state = initial.copy()
        self.time, self.steps = 0., 0
        self.euler_stages = 0
        self.source_work = self.ohmic_heat = self.outward_energy = 0.
        self.numerical_interface_heat = 0.
        self.ct_dissipation_heat = 0.  # poynting mode: estimate of sum E_diss*J (high-order fluxes)
        self.limited_euler_stages, self.minimum_high_order_weight = 0, 1.
        self.limiter_cell_counts = np.zeros((mesh.nr, mesh.nz), dtype=np.int64)
        self.limiter_troubled_counts = np.zeros_like(self.limiter_cell_counts)
        self.maximum_limiter_cell_fraction = 0.
        self.limiter_iterations = self.maximum_limiter_iterations = 0
        self.dual_energy_change = self.dual_energy_absolute_change = 0.
        self.dual_energy_cold_cell_steps = self.dual_energy_disagreements = 0
        self.dual_energy_max_cold_ratio = 0.
        self.dual_energy_discrepancy = 0.  # pressure mode: sum over cold cells of (E-K-M) - entropy thermal
        self.step_retries = 0
        self.floor_mass = self.floor_energy = 0.
        self.outflow_mass = 0.
        self.open_line_energy = 0.
        self.open_line_cell_steps = 0
        self.floor_cell_steps = 0
        self.hall_work = 0.             # time integral of the discrete sum J.E_Hall (zero to roundoff by construction)
        self.toroidal_boundary_work = 0.  # B_theta Poynting energy delivered through the gas boundaries (linked current)
        self.hall_limited_steps = self.hall_cutoff_cell_steps = 0
        self._hall_last = (False, 0)
        if initial.conserved.shape[-1] != (6 if p.toroidal else 4):
            raise ValueError("Toroidal runs need six conserved variables (with L and B_theta/r), others four")
        if p.toroidal_boundary is not None and not mesh.annular:
            raise ValueError("Only an annular chamber can be threaded by a current: the axis has r B_theta = 0")
        q0 = self._primitive(self.state, 0.)
        if p.dual_energy and self.state.entropy is None:
            self.state.entropy = entropy_density(q0[..., 0], q0[..., 3], p.gamma)
        elif not p.dual_energy:
            self.state.entropy = None
        u = self.mapping.interior(initial.induced_psi)
        ku = self.mapping.k@u
        if np.linalg.norm(ku[self.mapping.vacuum]) > 1e-10*max(np.linalg.norm(ku), 1e-30):
            raise ValueError("Initial induced field violates the vacuum constraint")
        if np.any(initial.induced_psi[0]):
            raise ValueError("Axis induced flux must vanish")
        if not self.mapping.open_boundary and (np.any(initial.induced_psi[-1]) or np.any(initial.induced_psi[:, [0, -1]])):
            raise ValueError("Fixed outer induced flux must vanish")
        self.initial = self._integrals(self.state, 0.)

    def _linked(self, time):
        """r B_theta of the vacuum beyond the gas boundaries at `time` (0 without a threading current)."""
        p = self.parameters
        return 0. if p.toroidal_boundary is None else float(p.toroidal_boundary(time))

    def total_flux(self, state=None, time=None):
        state = self.state if state is None else state
        time = self.time if time is None else time
        return state.induced_psi+self.sources.at(time).psi

    def _primitive(self, state, time):
        p = self.parameters
        flux = self.mapping.fluid_flux(self.total_flux(state, time))
        override = None
        if p.dual_energy and p.dual_energy_mode == "pressure" and state.entropy is not None:
            cold = self._cold(state.conserved, flux, state.entropy)
            override = (cold, state.entropy*state.conserved[..., 0]**(p.gamma-1))
        return primitive(self.mesh, MHDState(state.conserved, flux), p.gamma,
                         energy_quadrature=p.energy_quadrature, pressure_override=override)

    def _cold(self, conserved, flux, entropy):
        p = self.parameters
        return cold_cells(self.mesh, conserved, flux, p.gamma, p.dual_energy_threshold, p.energy_quadrature,
                          p.dual_energy_shock_threshold, p.dual_energy_shock_min_fraction, entropy,
                          p.dual_energy_criterion)[0]

    def _integrals(self, state, time):
        q = self._primitive(state, time)
        vol, m = self.mesh.volumes, self.mapping
        flux = m.fluid_flux(self.total_flux(state, time))
        magnetic_cell = float(np.sum(selected_energy(self.mesh, flux, self.parameters.energy_quadrature)*vol))
        energy = float(np.sum(state.conserved[..., 3]*vol))
        u = m.interior(state.induced_psi)
        induced_energy = float(np.pi/MU0*(u@(m.k@u)))
        result = {"mass_kg": float(np.sum(q[..., 0]*vol)), "fluid_total_energy_j": energy,
                  "gas_energy_j": energy-magnetic_cell, "fluid_magnetic_energy_j": magnetic_cell,
                  "induced_field_energy_j": induced_energy, "response_energy_j": energy-magnetic_cell+induced_energy,
                  "internal_energy_j": float(np.sum(q[..., 3]*vol/(self.parameters.gamma-1))),
                  "minimum_density_kg_m3": float(q[..., 0].min()), "minimum_pressure_pa": float(q[..., 3].min()),
                  "maximum_speed_m_s": float(np.hypot(q[..., 1], q[..., 2]).max())}
        if self.parameters.toroidal:   # gas_energy_j includes the B_theta energy: no B_theta exists outside the gas
            kinetic_t, magnetic_t = toroidal_energies(self.mesh, state.conserved)
            _, jz, jr = toroidal_ohm.currents(self.mesh, q[..., 7], self._linked(time))
            r = np.sqrt(toroidal_ohm.geometry(self.mesh).r2)
            result.update(angular_momentum_kg_m2_s=float(np.sum(state.conserved[..., 4]*vol)),
                          toroidal_flux_wb=float(np.sum(q[..., 7]*vol)/(2*np.pi)),
                          toroidal_kinetic_energy_j=float(np.sum(kinetic_t*vol)),
                          toroidal_magnetic_energy_j=float(np.sum(magnetic_t*vol)),
                          maximum_toroidal_field_t=float(np.abs(r*q[..., 7]).max()),
                          maximum_toroidal_speed_m_s=float(np.abs(r*q[..., 6]).max()),
                          maximum_axial_current_a_m2=float(np.abs(jz).max()),
                          maximum_radial_current_a_m2=float(np.abs(jr).max()))
        return result

    def timestep_limit(self):
        q = self._primitive(self.state, self.time)
        mesh, p = self.mesh, self.parameters
        if p.toroidal:
            r = np.sqrt(toroidal_ohm.geometry(mesh).r2)
            waves = np.max((abs(q[..., 1])+toroidal_fast_speed(q, r, 0, p.gamma))/mesh.dr+
                           (abs(q[..., 2])+toroidal_fast_speed(q, r, 1, p.gamma))/mesh.dz)
        else:
            waves = np.max((abs(q[..., 1])+fast_speed(q, 0, p.gamma))/mesh.dr+
                           (abs(q[..., 2])+fast_speed(q, 1, p.gamma))/mesh.dz)
        eta = interface_resistivity(self.mapping, q, p.gamma, p.interface_speed_cap) if p.interface_diffusion else np.zeros(1)
        physical = float(np.max(self._node_resistivity(q, self._current_magnitude(self.state, q, time=self.time))))
        diffusion = 2*(physical+float(eta.max()))/MU0*(1/mesh.dr**2+1/mesh.dz**2)
        hall = 0.
        if p.hall:
            # Whistler: omega = c|B|/mu0 k k_par <= 4 D (1/dr^2+1/dz^2) on the grid, D = max c|B|/mu0 over cells;
            # undamped, so RK3 needs omega dt <= sqrt(3). The rate enters the CFL sum like the other terms.
            coefficient, weight = self._hall_coefficient(q)
            b = np.sqrt(q[..., 4]**2+q[..., 5]**2+(r*q[..., 7])**2)
            omega = 4*float(np.max(coefficient*b))/MU0*(1/mesh.dr**2+1/mesh.dz**2)
            hall = p.cfl*omega/(p.hall_whistler_fraction*np.sqrt(3))
            self._hall_last = (bool(hall > waves+diffusion), int(np.sum(weight < .5)))
        limit = float(p.cfl/(waves+diffusion+hall))
        return limit if p.max_timestep is None else min(limit, p.max_timestep)

    def _hall_coefficient(self, q):
        p = self.parameters
        return toroidal_ohm.hall_coefficient(q[..., 0], p.hall_ion_mass or p.spitzer_ion_mass, p.hall_cutoff_density)

    def _current_magnitude(self, state, q, current=None, time=None):
        """Node current for the Chodura drift: |J_theta|, or |J| including J_r, J_z in 2.5D."""
        if current is None:
            current = self._node_current(state)[1]
        if not self.parameters.toroidal or self.parameters.chodura_coefficient == 0:
            return current
        m, mesh = self.mapping, self.mesh
        _, jz, jr = toroidal_ohm.currents(mesh, q[..., 7], self._linked(self.time if time is None else time))
        zn, rn = toroidal_ohm.node_currents(jz, jr)
        data = np.zeros(m.field.shape)
        data[m.r_offset:m.r_offset+mesh.nr+1, m.z_offset:m.z_offset+mesh.nz+1] = zn**2+rn**2
        return np.sqrt(current**2+m.grid_vector(data))

    def _fluid_node_resistivity(self, node_eta):
        """Physical eta on the fluid node grid (nr+1, nz+1); the axis row (not a mapping unknown) copies row 1."""
        mesh, m = self.mesh, self.mapping
        if np.ndim(node_eta) == 0:
            return np.full((mesh.nr+1, mesh.nz+1), float(node_eta))
        grid = np.array(m.fluid_flux(m.full(node_eta)))
        if not mesh.annular:
            grid[0] = grid[1]
        return grid

    def _toroidal_electric(self, state, q, current, node_eta, total, time):
        """Resistive and Hall E_z/E_r for B_theta/r, the Hall E_theta at fluid nodes, and their powers."""
        mesh, m, p = self.mesh, self.mapping, self.parameters
        _, jz, jr = toroidal_ohm.currents(mesh, q[..., 7], self._linked(time))
        ez, er, heat = toroidal_ohm.ohmic(mesh, self._fluid_node_resistivity(node_eta), jz, jr)
        etheta, hall_power = np.zeros((mesh.nr+1, mesh.nz+1)), 0.
        if p.hall:
            jtheta = np.array(m.fluid_flux(m.full(current)))
            if not mesh.annular:
                jtheta[0] = 0.   # J_theta vanishes on the axis
            brf, bzf = mesh.face_fields(total)
            etheta, hz, hr = toroidal_ohm.hall_fields(mesh, q[..., 7], jz, jr, jtheta, brf, bzf,
                                                     self._hall_coefficient(q)[0])
            if not mesh.annular:
                etheta[0] = 0.
            ez, er = ez+hz, er+hr
            hall_power = toroidal_ohm.work(mesh, jtheta, etheta, jz, hz, jr, hr)
        return (ez, er, ez[0]), etheta, heat, hall_power, (jz, jr)

    def _node_resistivity(self, q, current=None):
        """Physical resistivity at mapping nodes.

        Base: uniform scalar, or NRL perpendicular Spitzer 1.03e-4 lnL Te^-1.5 ohm m with
        Te = p/(2 n e) (Z=1, Te=Ti), cell values averaged to nodes. Optional Chodura anomalous
        term (theta-pinch formation standard): c m_e w_pi/(n e^2) [1-exp(-v_d/(f c_s))], with
        node density/sound speed averaged from cells and v_d = |J|/(n e) at the node.
        """
        p = self.parameters
        from scipy.constants import e as QE, m_e as ME, epsilon_0 as EPS0
        to_nodes = lambda cell: np.asarray(self.mapping.heat_to_cells.T@cell.ravel())
        if p.resistivity_model == "uniform":
            base = p.resistivity
        else:
            te = np.maximum(q[..., 3]/(2*q[..., 0]/p.spitzer_ion_mass*QE), p.spitzer_min_te_ev)
            base = to_nodes(1.03e-4*p.spitzer_coulomb_log*te**-1.5)
        if p.gap_resistivity > 0:
            rho = to_nodes(q[..., 0])
            base = base+p.gap_resistivity/(1+(rho/p.gap_density)**4)
        if p.chodura_coefficient == 0:
            return base
        if current is None:
            current = self._current_magnitude(self.state, q)
        n = to_nodes(q[..., 0]/p.spitzer_ion_mass)
        cs = to_nodes(np.sqrt(p.gamma*q[..., 3]/q[..., 0]))
        safe = np.where(n > 0, n, 1.)
        drift = np.abs(current)/(safe*QE)
        w_pi = np.sqrt(safe*QE**2/(EPS0*p.spitzer_ion_mass))
        anomalous = p.chodura_coefficient*ME*w_pi/(safe*QE**2)*(1-np.exp(-drift/(p.chodura_f*np.where(cs > 0, cs, 1.))))
        return base+np.where(n > 0, anomalous, 0.)

    def _node_current(self, state):
        m = self.mapping
        ku = m.k@m.interior(state.induced_psi)
        return ku, np.divide(ku, MU0*m.node_area, out=np.zeros(m.node_count), where=m.active)

    def _entropy_rate(self, state, time):
        """Advected rho*A with physical Ohmic heating only (numerical interface heat excluded).
        2.5D adds the paired eta (J_z^2 + J_r^2) face heat; the Hall field does no work."""
        m, p = self.mapping, self.parameters
        q = self._primitive(state, time)
        current = self._node_current(state)[1]
        node_eta = self._node_resistivity(q, self._current_magnitude(state, q, current, time))
        joule = m.fluid_flux(m.full(node_eta*current**2))
        heating = .25*(joule[:-1, :-1]+joule[1:, :-1]+joule[:-1, 1:]+joule[1:, 1:])
        if p.toroidal:
            geo = toroidal_ohm.geometry(self.mesh)
            _, jz, jr = toroidal_ohm.currents(self.mesh, q[..., 7], self._linked(time))
            ez, er, _ = toroidal_ohm.ohmic(self.mesh, self._fluid_node_resistivity(node_eta), jz, jr)
            heating = heating+toroidal_ohm.face_heat_to_cells(self.mesh, ez*jz*geo.radial, er*jr*geo.axial)
        return entropy_rate(self.mesh, q, state.entropy, heating, p.gamma, p.open_ends)

    def _fluxes(self, state, time, spatial_order=None):
        m, p = self.mapping, self.parameters
        q = self._primitive(state, time)
        ku, current = self._node_current(state)
        node_eta = self._node_resistivity(q, self._current_magnitude(state, q, current, time))
        electric = m.fluid_flux(m.full(node_eta*current))
        full_total = self.total_flux(state, time)
        total = m.fluid_flux(full_total)
        toroidal_electric, toroidal_heat, hall_power = None, 0., 0.
        if p.toroidal:
            toroidal_electric, etheta, toroidal_heat, hall_power, _ = self._toroidal_electric(state, q, current,
                                                                                          node_eta, total, time)
            if p.hall:
                electric = electric+etheta
        numerical_eta = interface_resistivity(m, q, p.gamma, p.interface_speed_cap) if p.interface_diffusion else np.zeros(m.node_count)
        numerical_electric = m.fluid_flux(m.full(numerical_eta*current))
        boundary = variational_trace(m, full_total, ku) if p.boundary_poynting == "variational" else None
        if boundary is not None:
            traces = boundary.face_fields(self.mesh)
        elif p.analytic_source_traces:
            traces = tangential_field(m, state.induced_psi, self.sources.at(time))
        else:
            traces = tangential_field(m, full_total)
        fluxes = flux_components(self.mesh, q, total, p.gamma, electric, traces,
                                  numerical_electric, p.spatial_order if spatial_order is None else spatial_order,
                                  p.local_trace_bounds, p.open_ends, p.magnetic_energy_flux, toroidal_electric,
                                  self._linked(time) if p.toroidal else 0.)
        fluxes.boundary_poynting = boundary
        source_derivative = m.grid_vector(self.sources.at(time).dpsi_dt)
        source_power = -float((2*np.pi/MU0)*ku[m.active]@source_derivative[m.active])
        ohmic_power = float(2*np.pi*np.sum(node_eta*current**2*m.node_radius*m.node_area))+toroidal_heat
        fluxes.hall_power = hall_power
        numerical_power = float(2*np.pi*np.sum(numerical_eta*current**2*m.node_radius*m.node_area))
        fluxes.ct_power = 0.
        if fluxes.dissipative_electric is not None:
            weight = m.fluid_flux(m.full(current*m.node_radius*m.node_area))
            fluxes.ct_power = float(2*np.pi*np.sum(fluxes.dissipative_electric*weight))
        return fluxes, (source_power, ohmic_power, numerical_power)

    def _forward_euler(self, state, time, total, dt):
        """Select shared face/corner fluxes before assembly; never repair energy."""
        high, powers = self._fluxes(state, time)
        p = self.parameters
        mode = p.flux_limiter if p.positivity_flux_limiting and p.spatial_order == 2 else "off"
        exempt = self._cold(state.conserved, total, state.entropy) if p.dual_energy else None
        rhs, report = admissible_update(self.mesh, state.conserved, total, dt, high,
                                        lambda: self._fluxes(state, time, spatial_order=1)[0], mode=mode,
                                        energy_quadrature=p.energy_quadrature, exempt=exempt)
        entropy = self._entropy_rate(state, time) if p.dual_energy else None
        outflow_mass = -float(np.sum(rhs[0][..., 0]*self.mesh.volumes))  # interior fluxes telescope
        # A linked vacuum B_theta adds the boundary Poynting inflow as a 4th rhs entry (actually selected fluxes).
        linked_inflow = rhs[3] if len(rhs) > 3 else 0.
        # Current, physical/interface eta and analytic source power are independent
        # of reconstruction order at this state; only outward flux is selected.
        return (*rhs[:3], *powers, outflow_mass, entropy, high.ct_power, high.hall_power, linked_inflow), report

    def _new_state(self, conserved, fluid_total_flux, time, entropy=None, final=False):
        m, mesh, p = self.mapping, self.mesh, self.parameters
        data = np.zeros(m.field.shape)
        external = self.sources.at(time).psi
        data[m.r_offset:m.r_offset+mesh.nr+1, m.z_offset:m.z_offset+mesh.nz+1] = fluid_total_flux-m.fluid_flux(external)
        if p.dual_energy:
            if np.any(conserved[..., 0] <= 0) or not np.all(np.isfinite(entropy)) or np.any(entropy <= 0):
                raise NonphysicalState("Dual-energy stage has nonpositive density or entropy; no repair")
            chosen = select(mesh, conserved, fluid_total_flux, entropy, p.gamma, p.dual_energy_threshold,
                            p.energy_quadrature, p.dual_energy_shock_threshold, p.dual_energy_shock_min_fraction,
                            p.dual_energy_criterion)
            entropy = chosen.entropy
            if p.dual_energy_mode == "reset":
                conserved = chosen.conserved
        result = ConservativeState(conserved, m.harmonic_extension(data), entropy)
        self._primitive(result, time)
        if p.dual_energy and final:
            self._pending_dual = chosen
        return result

    def step(self, dt=None):
        limit = self.timestep_limit()
        dt = limit if dt is None else dt
        if not np.isfinite(dt) or dt <= 0 or dt > limit*(1+1e-12):
            raise ValueError(f"Time step must be positive and <= {limit:.6g} s")
        for attempt in range(self.parameters.max_step_halvings+1):
            try:
                return self._attempt(dt)
            except NonphysicalState:
                if attempt == self.parameters.max_step_halvings:
                    raise
                dt *= .5
                self.step_retries += 1

    def _floored(self, state, time):
        """Apply the declared minimum-density closure (Lie split after the MHD step) and revalidate."""
        p = self.parameters
        flux = self.mapping.fluid_flux(self.total_flux(state, time))
        if p.dual_energy:
            cold = self._cold(state.conserved, flux, state.entropy)
        else:
            cold = np.zeros(state.conserved.shape[:2], dtype=bool)
        floor = apply_floor(self.mesh, state.conserved, state.entropy, self._cell_magnetic(state, flux),
                            cold, p.minimum_density, p.floor_temperature_ev, p.floor_ion_mass, p.gamma)
        result = ConservativeState(floor.conserved, state.induced_psi, floor.entropy)
        self._primitive(result, time)
        return result, floor

    def _cell_magnetic(self, state, flux):
        """Cell magnetic energy density: the selected poloidal quadrature plus, in 2.5D, B_theta^2/(2 mu0)."""
        magnetic = selected_energy(self.mesh, flux, self.parameters.energy_quadrature)
        if self.parameters.toroidal:
            magnetic = magnetic+toroidal_energies(self.mesh, state.conserved)[1]
        return magnetic

    def _open_line_loss(self, state, time, dt):
        """Declared flux-limited parallel conduction loss on open field lines (psi >= 0), Lie split."""
        p = self.parameters
        flux = self.mapping.fluid_flux(self.total_flux(state, time))
        if p.dual_energy:
            cold = self._cold(state.conserved, flux, state.entropy)
        else:
            cold = np.zeros(state.conserved.shape[:2], dtype=bool)
        psi_cells = .25*(flux[:-1, :-1]+flux[1:, :-1]+flux[:-1, 1:]+flux[1:, 1:])
        loss = open_line_loss(self.mesh, state.conserved, state.entropy, self._cell_magnetic(state, flux),
                              cold, psi_cells, dt, p.spitzer_ion_mass, p.gamma, p.open_line_flux_limit,
                              p.open_line_wall_temperature_ev)
        result = ConservativeState(loss.conserved, state.induced_psi, loss.entropy)
        self._primitive(result, time)
        return result, loss

    def _attempt(self, dt):
        """One SSP step from the accepted state; nothing is committed unless every stage validates."""
        old = self.state
        total = self.mapping.fluid_flux(self.total_flux())
        s0 = old.entropy
        # Entropy uses one forward-Euler first-order step from the start state (positive under CFL);
        # every stage sees the same advanced rho*A. Warm cells resynchronize it from total energy.
        combine = (lambda *terms: s0+dt*a[7]) if s0 is not None else (lambda *terms: None)
        a, report_a = self._forward_euler(old, self.time, total, dt)
        predictor = self._new_state(old.conserved+dt*a[0], total+dt*a[1], self.time+dt, combine((1, a)))
        b, report_b = self._forward_euler(predictor, self.time+dt, total+dt*a[1], dt)
        if self.parameters.time_order == 2:
            final = self._new_state(old.conserved+.5*dt*(a[0]+b[0]), total+.5*dt*(a[1]+b[1]), self.time+dt,
                                    combine((.5, a), (.5, b)), final=True)
            weighted = ((.5, a), (.5, b))
            reports = (report_a, report_b)
        else:
            half_flux = total+.25*dt*(a[1]+b[1])
            half = self._new_state(old.conserved+.25*dt*(a[0]+b[0]), half_flux, self.time+.5*dt,
                                   combine((.25, a), (.25, b)))
            c, report_c = self._forward_euler(half, self.time+.5*dt, half_flux, dt)
            final = self._new_state(old.conserved+dt*((a[0]+b[0])/6+2*c[0]/3),
                                    total+dt*((a[1]+b[1])/6+2*c[1]/3), self.time+dt,
                                    combine((1/6, a), (1/6, b), (2/3, c)), final=True)
            weighted = ((1/6, a), (1/6, b), (2/3, c))
            reports = (report_a, report_b, report_c)
        floor = loss = None
        if self.parameters.minimum_density is not None:
            final, floor = self._floored(final, self.time+dt)
        if self.parameters.open_line_flux_limit > 0:
            final, loss = self._open_line_loss(final, self.time+dt, dt)
        # Every Euler proposal and exact vacuum completion validates atomically.
        self.state = final
        if loss is not None:
            self.open_line_energy += loss.energy_j
            self.open_line_cell_steps += loss.cells
        if floor is not None:
            self.floor_mass += floor.mass_kg
            self.floor_energy += floor.energy_j
            self.floor_cell_steps += floor.cells
        if self.parameters.dual_energy:
            chosen, p = self._pending_dual, self.parameters
            if p.dual_energy_mode == "reset":
                self.dual_energy_change += chosen.energy_change_j
                self.dual_energy_absolute_change += chosen.absolute_change_j
            else:  # total energy untouched; current (E-K-M) minus entropy thermal energy over cold cells
                self.dual_energy_discrepancy = -chosen.energy_change_j
            self.dual_energy_cold_cell_steps += int(chosen.cold.sum())
            self.dual_energy_disagreements += chosen.disagreements
            self.dual_energy_max_cold_ratio = max(self.dual_energy_max_cold_ratio, chosen.max_cold_ratio)
        self.time += dt; self.steps += 1
        self.euler_stages += len(reports)
        self.outward_energy += dt*sum(w*rates[2] for w, rates in weighted)
        self.source_work += dt*sum(w*rates[3] for w, rates in weighted)
        self.ohmic_heat += dt*sum(w*rates[4] for w, rates in weighted)
        self.numerical_interface_heat += dt*sum(w*rates[5] for w, rates in weighted)
        self.ct_dissipation_heat += dt*sum(w*rates[8] for w, rates in weighted)
        if self.parameters.toroidal_boundary is not None:
            self.toroidal_boundary_work += dt*sum(w*rates[10] for w, rates in weighted)
        if self.parameters.hall:
            self.hall_work += dt*sum(w*rates[9] for w, rates in weighted)
            self.hall_limited_steps += int(self._hall_last[0])
            self.hall_cutoff_cell_steps += self._hall_last[1]
        self.outflow_mass += dt*sum(w*rates[6] for w, rates in weighted)
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
        # Reset-mode dual-energy corrections (zero in pressure mode) are a declared, counted source.
        response_residual = (result["response_energy_j"]-self.initial["response_energy_j"]-self.source_work-self.floor_energy
                             +self.open_line_energy-self.dual_energy_change)
        if self.parameters.toroidal_boundary is not None:   # work of the threading current, counted separately
            response_residual -= self.toroidal_boundary_work
        fluid_residual = (result["fluid_total_energy_j"]-self.initial["fluid_total_energy_j"]+self.outward_energy-self.floor_energy
                          +self.open_line_energy-self.dual_energy_change)
        ku = self.mapping.k@self.mapping.interior(self.state.induced_psi)
        result.update(time_s=self.time, steps=self.steps, euler_stages=self.euler_stages,
                      source_work_j=self.source_work, ohmic_heat_j=self.ohmic_heat,
                      numerical_interface_heat_j=self.numerical_interface_heat,
                      ct_dissipation_heat_j=self.ct_dissipation_heat,
                      limited_euler_stages=self.limited_euler_stages, minimum_high_order_weight=self.minimum_high_order_weight,
                      affected_euler_cell_stages=int(self.limiter_cell_counts.sum()),
                      troubled_euler_cell_stages=int(self.limiter_troubled_counts.sum()),
                      maximum_limiter_cell_fraction=self.maximum_limiter_cell_fraction,
                      limiter_iterations=self.limiter_iterations, maximum_limiter_iterations=self.maximum_limiter_iterations,
                      dual_energy_change_j=self.dual_energy_change,
                      dual_energy_absolute_change_j=self.dual_energy_absolute_change,
                      dual_energy_cold_cell_steps=self.dual_energy_cold_cell_steps,
                      dual_energy_disagreement_cell_steps=self.dual_energy_disagreements,
                      dual_energy_max_cold_ratio=self.dual_energy_max_cold_ratio,
                      dual_energy_cold_discrepancy_j=self.dual_energy_discrepancy, step_retries=self.step_retries,
                      outward_fluid_energy_j=self.outward_energy, fluid_energy_balance_residual_j=fluid_residual,
                      response_energy_balance_residual_j=response_residual,
                      response_residual_over_source_work=response_residual/abs(self.source_work) if self.source_work else None,
                      response_energy_balance_relative=response_residual/max(abs(self.initial["response_energy_j"]), abs(self.source_work), 1e-30),
                      mass_relative_change=(result["mass_kg"]-self.initial["mass_kg"]-self.floor_mass+self.outflow_mass)/self.initial["mass_kg"],
                      outflow_mass_kg=self.outflow_mass, open_line_loss_j=self.open_line_energy,
                      hall_work_j=self.hall_work, hall_limited_steps=self.hall_limited_steps,
                      toroidal_boundary_work_j=self.toroidal_boundary_work,
                      hall_cutoff_cell_steps=self.hall_cutoff_cell_steps,
                      open_line_cell_steps=self.open_line_cell_steps,
                      floor_mass_kg=self.floor_mass, floor_energy_j=self.floor_energy, floor_cell_steps=self.floor_cell_steps,
                      vacuum_relative_residual=float(np.linalg.norm(ku[self.mapping.vacuum])/max(np.linalg.norm(ku), 1e-30)),
                      max_face_divergence_t_m=float(abs(self.mesh.divergence(self.mapping.fluid_flux(self.total_flux()))).max()))
        return result
