"""Coarse NUCTE-III field-reversed theta-pinch formation on the fixed-chamber conservative MHD path.

Defaults are the benchmark configuration of docs/NUCTE_BENCHMARK.md (case-1 coil); case 6 is
end_radius=0.1267 filaments=12 bias_t=-0.09."""
import json
import sys
import time
from dataclasses import dataclass, asdict, replace
from pathlib import Path
import numpy as np
from scipy.constants import k as KB, e as QE, atomic_mass as AMU
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.geometry import Grid
from plasma_bench.coils import Coil
from plasma_bench.waveforms import Callback
from plasma_bench.electromagnetics import VacuumField, loop_field
from plasma_bench.mhd_state import from_primitive
from plasma_bench.conservative_coupled import ConservativeVacuumMHD, ConservativeParameters, ConservativeState
from plasma_bench.coil_circuit import Strap, section_coil, Module, Crowbar, simulate, CircuitStepper
from scipy.constants import mu_0 as MU0

D_MASS = 2.01410*AMU


@dataclass(frozen=True)
class Nucte:
    """Published NUCTE-III geometry/fill (Takahashi 2000; Kuroyanagi 2014) plus explicit assumptions."""
    wall_radius: float = .115         # quartz tube, 23 cm ID, 2.0 m long
    half_length: float = 1.0
    coil_radius: float = .17          # 18 single-turn straps, 5 cm wide
    end_radius: float = .162          # 5 straps each end; (0.17/0.162)^2 = 1.10, case 1 mirror ratio
    taper_start: float = 0.           # >0: graded coil, radius linear from coil_radius at |z|=taper_start to
                                      # end_radius at the coil end (NUCTE-III stepped 36..26 cm elements)
    strap_width: float = .05
    centre_straps: int = 18
    end_straps: int = 5
    filaments: int = 4                # per strap, equal share of the strap current
    fill_pa: float = 1.33             # D2
    gas_k: float = 300.
    pi_temperature_ev: float = 2.     # ASSUMED complete pre-ionization, Te = Ti
    pi_half_length: float = .7        # >0: ionized only for |z| < this (theta-PI under the coil); beyond, floor-level plasma
    pi_edge_width: float = .03        # tanh transition width of the pre-ionized region [m]
    bias_t: float = -.06              # ASSUMED; published maximum is 0.09 T
    peak_t: float = .4                # ASSUMED vacuum centre field at end of rise
    quarter_period_s: float = 4e-6    # published rise time
    crowbar_decay_s: float = 100e-6   # ASSUMED L/R after crowbar
    resistivity: float = 1e-5         # ASSUMED uniform, Ohm m (Spitzer at ~13 eV); used when resistivity_model="uniform"
    resistivity_model: str = "spitzer"  # or "uniform": NRL eta_perp of local Te (Te=Ti, lnL=10, Te>=1 eV for eta)
    chodura: float = .1               # Chodura anomalous-resistivity coefficient (standard ~0.1; 0 disables)
    nr: int = 16
    nz: int = 128
    exterior_cells: int = 6           # fixed induced-flux boundary just inside the coil
    end_cells: int = 2                # vacuum field cells beyond each tube end before the fixed boundary
    end_s: float = 15e-6
    cfl: float = .2
    dual_energy: bool = True          # entropy pressure in cold/low-beta cells (pressure mode, energy exact)
    dual_energy_mode: str = "reset"   # cold-cell total energy reset to the entropy value (counted); or "pressure"
    dual_energy_criterion: str = "entropy"  # cold test on the entropy thermal energy; "total" (E-K-M) can make
                                            # hot cells unbounded reservoirs in pressure mode (dual_energy.md)
    max_step_halvings: int = 6        # timestep control: retry inadmissible steps with dt/2
    interface_diffusion: bool = True  # solver's declared numerical resistivity at the chamber boundary
    interface_speed_cap: float = 0.   # >0 caps the signal speed in that numerical resistivity [m/s]
    open_ends: bool = False           # outflow-only tube ends (open field-line end losses) instead of reflecting
    conduction_limit: float = 0.      # >0: flux-limited parallel conduction loss on open lines (f ~ 0.1)
    energy_quadrature: str = "staggered"  # with variational boundary Poynting flux and RK3: the compatible
    boundary_poynting: str = "variational"  # magnetic-energy path (docs/COMPATIBLE_ENERGY.md); closes the
    time_order: int = 3                     # response energy to ~1e-4 (historical: cell/riemann/2)
    magnetic_energy_flux: str = "ct"  # magnetic energy moved by the CT field's Poynting flux; historical "face"
    analytic_traces: bool = True      # wall traces: differenced induced flux + analytic coil field (needed when
                                      # coil filaments sit within two field cells of the wall, e.g. case 6)
    local_trace_bounds: bool = True   # reconstructed density/pressure bounded by neighbour minima
    floor_fraction: float = .01       # DECLARED gap closure: minimum density as a fraction of the fill (0 disables)
    drive: str = "ideal"              # "ideal": centre_field waveform; "circuit": capacitor banks through coil_circuit,
                                      # voltages set so the vacuum centre field is bias_t at t=0 and peaks at peak_t;
                                      # "coupled": same banks, but the plasma's image current loads them each step
    sections: int = 1                 # circuit: axial sections of the main coil (1 = common collector)
    main_modules: int = 1             # circuit: main modules per section
    main_capacitance: float = 67.5e-6  # circuit: total main-bank capacitance [F] (NUCTE-III fast bank)
    main_stray: float = 27e-9         # circuit: main-bank stray inductance seen by the coil, all modules in parallel [H]
    main_resistance: float = 0.       # circuit: main-bank series resistance, all modules in parallel [ohm]
    crowbar_inductance: float = 5e-9  # circuit: crowbar branch inductance, all crowbars in parallel [H]
    unidirectional_start: bool = True  # circuit: main start switches open at their first current zero (ignitron-like)
    section_skew_s: float = 0.        # circuit: sections centred at z > 0 fire this much later [s]
    bias_capacitance: float = 2000e-6  # circuit: bias bank (NUCTE-III, Kuroyanagi 2014)
    bias_inductance: float = 2.4e-6
    bias_resistance: float = 8e-3


def strap_layout(c):
    w = c.strap_width
    centre = [(c.coil_radius, (k-(c.centre_straps-1)/2)*w) for k in range(c.centre_straps)]
    edge = c.centre_straps*w/2
    straps = centre+[(c.end_radius, s*(edge+(k+.5)*w)) for s in (-1, 1) for k in range(c.end_straps)]
    if c.taper_start > 0:
        coil_end = edge+c.end_straps*w
        def radius(z):
            frac = np.clip((abs(z)-c.taper_start)/(coil_end-c.taper_start), 0, 1)
            return c.coil_radius+(c.end_radius-c.coil_radius)*frac
        straps = [(radius(z), z) for _, z in straps]
    return straps


def filaments(c, strap):
    r, z = strap
    return [(r, z+((j+.5)/c.filaments-.5)*c.strap_width) for j in range(c.filaments)]


def equal_flux_currents(c, inset=.003):
    """Parallel straps on one collector link equal flux; solve M I = 1, normalize to 1 T at centre."""
    fil = [filaments(c, s) for s in strap_layout(c)]
    n = len(fil)
    m = np.zeros((n, n))
    for a, fa in enumerate(fil):
        za = np.array([z for _, z in fa])
        ra = np.full(za.size, fa[0][0]-inset)
        for b, fb in enumerate(fil):
            m[a, b] = 2*np.pi*np.mean(sum(loop_field(ra, za, r, z)[0] for r, z in fb)/len(fb))
    current = np.linalg.solve(m, np.ones(n))
    b0 = sum(current[k]*sum(float(loop_field(0., 0., r, z)[2]) for r, z in fil[k])/len(fil[k]) for k in range(n))
    return fil, current/b0


def centre_field(c):
    tq, tc, amp = c.quarter_period_s, c.crowbar_decay_s, c.peak_t-c.bias_t

    def value(t):
        t = np.asarray(t, float)
        rise = c.bias_t+amp*np.sin(np.pi*np.clip(t, 0, tq)/(2*tq))
        return np.where(t <= tq, rise, c.peak_t*np.exp(-(t-tq)/tc))

    def rate(t):
        t = np.asarray(t, float)
        rise = amp*np.pi/(2*tq)*np.cos(np.pi*np.clip(t, 0, tq)/(2*tq))
        return np.where(t < 0, 0., np.where(t <= tq, rise, -c.peak_t/tc*np.exp(-(t-tq)/tc)))
    return Callback(value, rate, "bias, sinusoidal main rise to quarter period, crowbar L/R decay [T at centre]")


def coil_sections(c):
    """Main-coil straps assigned to `sections` contiguous axial groups (one collector each)."""
    layout = sorted(strap_layout(c), key=lambda s: s[1])
    per = int(np.ceil(len(layout)/c.sections))
    return section_coil([Strap(r, z, c.strap_width, section=k//per) for k, (r, z) in enumerate(layout)], c.filaments)


def circuit_drive(c, coil):
    """Bias and main banks (coil_circuit) on the main-coil sections. Linear superposition sets the charge
    voltages: bias so the vacuum centre field is bias_t at t=0, main so the total centre field peaks at
    peak_t (no-skew design values; skew then only shifts sections in time). Crowbars close at each main
    capacitor's voltage zero; their resistance gives the L/R decay crowbar_decay_s."""
    bz = np.array([f[2] for f in coil.field(0., 0.)], dtype=float)
    n_s = coil.sections
    n_mod = n_s*c.main_modules
    l_coil = float(1/np.linalg.inv(coil.inductance).sum())   # all sections in parallel
    r_cb = (l_coil+c.crowbar_inductance)/c.crowbar_decay_s
    centres = [np.mean([s.z for s in coil.straps if s.section == k]) for k in range(n_s)]

    def banks(v_main, v_bias, lead):
        mods = [Module(c.bias_capacitance/n_s, v_bias, c.bias_inductance*n_s, c.bias_resistance*n_s, section=k,
                       fire_time=-lead, name=f"bias{k}") for k in range(n_s)]
        mods += [Module(c.main_capacitance/n_mod, v_main, c.main_stray*n_mod, c.main_resistance*n_mod, section=k,
                        fire_time=c.section_skew_s if centres[k] > 0 else 0.,
                        crowbar=Crowbar(c.crowbar_inductance*n_mod, r_cb*n_mod), unidirectional=c.unidirectional_start,
                        name=f"main{k}.{j}")
                 for k in range(n_s) for j in range(c.main_modules)]
        return mods

    times = lambda lead: np.concatenate((np.linspace(-lead, 0, 401)[:-1], np.arange(0, c.end_s+1e-9, 2e-9)))
    centre = lambda res: res.section_current@bz
    # Bias alone at 1 V: lead time to its current peak, then the voltage for bias_t at t = 0.
    probe = simulate(coil.inductance, [m for m in banks(0., 1., 0.) if m.name.startswith("bias")], 400e-6, samples=8001)
    lead = float(probe.time[np.argmax(np.abs(probe.section_current).sum(axis=1))])
    v_bias = c.bias_t/float(np.interp(lead, probe.time, centre(probe)))
    v_main = (c.peak_t-c.bias_t)/float(centre(simulate(coil.inductance, [m for m in banks(1., 0., 0.) if m.name.startswith("main")],
                                                          2*c.quarter_period_s, samples=2001)).max())
    for _ in range(3):
        res = simulate(coil.inductance, banks(v_main, v_bias, lead), c.end_s, samples=times(lead))
        peak = float(centre(res)[res.time >= 0].max())
        v_main *= (c.peak_t-np.interp(0., res.time, centre(res)))/(peak-np.interp(0., res.time, centre(res)))
    res = simulate(coil.inductance, banks(v_main, v_bias, lead), c.end_s, samples=times(lead))
    mods = res.modules
    main = [k for k, m in enumerate(mods) if m.name.startswith("main")]
    summary = dict(section_inductance_nh=(coil.inductance*1e9).tolist(), coil_parallel_nh=l_coil*1e9,
                   centre_field_t_per_ka=(bz*1e3).tolist(), main_voltage_v=v_main, bias_voltage_v=v_bias,
                   bias_lead_s=lead, crowbar_resistance_per_module_ohm=r_cb*n_mod,
                   crowbar_times_s=[res.crowbar_times[k] for k in main],
                   main_peak_current_per_module_a=[float(np.abs(res.module_current[:, k]).max()) for k in main],
                   main_reversal_v=[float(res.capacitor_voltage[:, k].min()) for k in main],
                   main_energy_j=float(sum(.5*mods[k].capacitance*mods[k].voltage**2 for k in main)),
                   centre_field_t_at_0=float(np.interp(0., res.time, centre(res))),
                   centre_field_peak_t=float(centre(res)[res.time >= 0].max()),
                   energy_residual_relative=float(np.abs(res.energy_residual(coil.inductance)).max()
                                                  / sum(.5*m.capacitance*m.voltage**2 for m in mods)))
    return res, summary


def build(c):
    fluid = MHDMesh(c.wall_radius, c.half_length, c.nr, c.nz)
    grid = Grid(c.wall_radius+c.exterior_cells*fluid.dr, c.half_length+c.end_cells*fluid.dz,
                c.nr+c.exterior_cells+1, c.nz+2*c.end_cells+1)
    stepper = None
    if c.drive in ("circuit", "coupled"):
        coil = coil_sections(c)
        circuit, drive = circuit_drive(c, coil)
        if c.drive == "coupled":
            stepper = CircuitStepper(coil.inductance, circuit.modules)
            stepper.run_to(0., 20e-9)
            coils = coil.sources(stepper.waveforms)
        else:
            coils = coil.sources(circuit.section_waveforms())
        current = np.zeros(1)
    else:
        fil, current = equal_flux_currents(c)
        wave = centre_field(c)
        coils = tuple(Coil(r, z, wave, turns=float(i)/c.filaments, name=f"strap{k}")
                      for k, (strap, i) in enumerate(zip(fil, current)) for r, z in strap)
        drive = None
    sources = VacuumField(grid, coils)
    n_d = 2*c.fill_pa/(KB*c.gas_k)
    rho, p = n_d*D_MASS, 2*n_d*c.pi_temperature_ev*QE
    if c.pi_half_length > 0:
        # ASSUMED theta-PI profile: full ionization under the coil, residual (floor-level) plasma beyond
        z = fluid.z[None, :]
        shape = .5*(1-np.tanh((abs(z)-c.pi_half_length)/c.pi_edge_width))
        residual = max(c.floor_fraction, 1e-3)
        fraction = residual+(1-residual)*shape
        rho, p = rho*fraction*np.ones((c.nr, 1)), p*fraction*np.ones((c.nr, 1))
    params = ConservativeParameters(resistivity=c.resistivity, resistivity_model=c.resistivity_model, chodura_coefficient=c.chodura,
                                    spitzer_ion_mass=D_MASS, vacuum_boundary="fixed", cfl=c.cfl,
                                    dual_energy=c.dual_energy, max_step_halvings=c.max_step_halvings,
                                    dual_energy_mode=c.dual_energy_mode, dual_energy_criterion=c.dual_energy_criterion,
                                    energy_quadrature=c.energy_quadrature, boundary_poynting=c.boundary_poynting,
                                    time_order=c.time_order,
                                    interface_diffusion=c.interface_diffusion,
                                    interface_speed_cap=c.interface_speed_cap if c.interface_speed_cap > 0 else None,
                                    open_ends=c.open_ends, open_line_flux_limit=c.conduction_limit,
                                    analytic_source_traces=c.analytic_traces, magnetic_energy_flux=c.magnetic_energy_flux,
                                    open_line_wall_temperature_ev=c.pi_temperature_ev,
                                    local_trace_bounds=c.local_trace_bounds,
                                    minimum_density=c.floor_fraction*2*c.fill_pa/(KB*c.gas_k)*D_MASS if c.floor_fraction > 0 else None,
                                    floor_temperature_ev=c.pi_temperature_ev, floor_ion_mass=D_MASS)
    # Initial plasma is current-free in the bias: induced flux zero, total flux = source flux.
    probe = ConservativeState(np.zeros((c.nr, c.nz, 4)), np.zeros((grid.nr, grid.nz)))
    total = sources.at(0.).psi[:c.nr+1, c.end_cells:c.end_cells+c.nz+1]
    state = from_primitive(fluid, total, rho, p, energy_quadrature=c.energy_quadrature)
    model = ConservativeVacuumMHD(fluid, sources, ConservativeState(state.conserved, probe.induced_psi), params)
    model.coil_half_length = (c.centre_straps/2+c.end_straps)*c.strap_width
    model.coupler = None
    if stepper is not None:
        model.coupler = coupler(model, coil, stepper, sources)
    return model, dict(fill_deuterons_m3=n_d, density_kg_m3=n_d*D_MASS, pressure_pa=2*n_d*c.pi_temperature_ev*QE,
                       strap_currents_a_per_t=current.tolist(),
                       centre_strap_share=float(current[:c.centre_straps].sum()/current.sum()) if c.drive == "ideal" else None,
                       circuit=drive)


def separatrix(model):
    """r_s(z) where total flux returns to zero from a reversed core; NaN where the axis is unreversed."""
    m = model.mesh
    psi = model.mapping.fluid_flux(model.total_flux())
    rf = m.rf
    rs = np.full(psi.shape[1], np.nan)
    for j in range(psi.shape[1]):
        col = psi[:, j]
        if col[1] < 0:
            k = np.nonzero(col[1:] >= 0)[0]
            if k.size:
                i = k[0]+1
                rs[j] = rf[i-1]+(rf[i]-rf[i-1])*(-col[i-1])/(col[i]-col[i-1])
            else:
                rs[j] = np.inf  # reversed out to the wall: not separated from the wall
    return m.zf, rs, psi


def excluded_flux_radius(model, psi):
    """Experimental definition: r_dphi^2 = r_w^2 - flux_loop/(pi B_wall), at every axial face."""
    m = model.mesh
    b_wall = (psi[-1]-psi[-2])/(m.rf[-1]-m.rf[-2])*2/(m.rf[-1]+m.rf[-2])
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.sqrt(np.maximum(m.radius**2-2*psi[-1]/b_wall, 0.))


def section_image_currents(model, coil):
    """Image current the fixed induced-flux boundary carries over each coil section's axial extent [A]."""
    grid = model.mapping.field
    u = model.state.induced_psi
    bz = (-4*u[-2]+u[-3])/(2*grid.dr)/grid.radius
    out = np.zeros(coil.sections)
    for k in range(coil.sections):
        straps = [st for st in coil.straps if st.section == k]
        lo, hi = min(st.z-st.width/2 for st in straps), max(st.z+st.width/2 for st in straps)
        inside = (grid.z >= lo-1e-12) & (grid.z <= hi+1e-12)
        out[k] = np.trapezoid(np.where(inside, bz, 0.), grid.z)/MU0
    return out


def coupler(model, coil, stepper, sources):
    """on_step hook: advance the bank circuit to the accepted solver time with the plasma as its load,
    then refresh the coil waveforms (the vacuum field cache is cleared).

    Load: the plasma current's flux linkage with each section, Phi_p,s = sum_j (2 pi/mu0) (K u)_j psi_s(j)
    (psi_s: the section's vacuum flux per ampere), enters as the equivalent image current
    I_img = -L^-1 Phi_p, so V_s = d(L I_coil + Phi_p)/dt with I_coil = I_vac + I_img. The power the banks
    deliver through it, L dI_vac/dt . I_img = -dI_vac/dt . Phi_p, is exactly the solver's source work."""
    m = model.mapping
    units = [m.grid_vector(np.asarray(group[1][0])) for group in sources._groups]
    l_inv = np.linalg.inv(coil.inductance)

    def update(mm):
        ku = mm._node_current(mm.state)[0]
        phi = np.array([2*np.pi/MU0*float(ku[m.active]@unit[m.active]) for unit in units])
        stepper.advance(mm.time, -l_inv@phi)
        sources._cache.clear()
    update.stepper = stepper
    return update


def image_current(model, coil_half_length):
    """Current the fixed induced-flux boundary (the coil as a flux conserver) carries over the coil length:
    (1/mu0) integral of the induced Bz just inside r_b. A circuit with stray inductance supplies less."""
    grid = model.mapping.field
    u = model.state.induced_psi
    bz = (-4*u[-2]+u[-3])/(2*grid.dr)/grid.radius
    inside = np.abs(grid.z) <= coil_half_length
    return float(np.trapezoid(np.where(inside, bz, 0.), grid.z)/MU0)


def observe(model):
    z, rs, psi = separatrix(model)
    q = model._primitive(model.state, model.time)
    finite = np.isfinite(rs)
    closed = finite & (rs < model.mesh.radius)
    mid = len(z)//2
    rdphi = excluded_flux_radius(model, psi)
    return dict(time_s=model.time, steps=model.steps, axis_bz_min_t=float((psi[1]/(model.mesh.rf[1]**2/2)).min()),
                separatrix_mid_m=float(rs[mid]) if np.isfinite(rs[mid]) else None,
                reversed_length_m=float(np.ptp(z[finite])) if finite.sum() > 1 else 0.,
                separated_length_m=float(np.ptp(z[closed])) if closed.sum() > 1 else 0.,
                excluded_radius_mid_m=float(rdphi[mid]),
                axis_density_m3=float(q[0, q.shape[1]//2, 0]/D_MASS),
                max_density_m3=float(q[..., 0].max()/D_MASS), min_density_m3=float(q[..., 0].min()/D_MASS),
                max_pressure_pa=float(q[..., 3].max()), max_speed_m_s=float(np.hypot(q[..., 1], q[..., 2]).max()),
                image_current_a=image_current(model, model.coil_half_length),
                circuit=model.coupler.stepper.snapshot() if model.coupler else None)


def main(out=Path("outputs/nucte_formation"), **overrides):
    c = replace(Nucte(), **overrides)
    out.mkdir(parents=True, exist_ok=True)
    model, setup = build(c)
    (out/"setup.json").write_text(json.dumps(dict(case=asdict(c), **setup), indent=1))
    history, snapshots, velocity, start = [observe(model)], {}, {}, time.time()
    marks = list(np.arange(1e-6, c.end_s+1e-12, 1e-6))
    failure = None
    while model.time < c.end_s:
        target = marks.pop(0) if marks else c.end_s
        try:
            model.advance(target, on_step=model.coupler)
        except Exception as error:  # record the accepted state and the failure; never repaired
            failure = dict(time_s=model.time, steps=model.steps, error=f"{type(error).__name__}: {error}")
            (out/"failure.json").write_text(json.dumps(failure, indent=1))
            print(json.dumps(failure), flush=True)
            break
        row = observe(model)
        row.update(wall_s=time.time()-start, **{k: v for k, v in model.diagnostics().items()
                                                if k in ("mass_relative_change", "response_energy_balance_relative",
                                                         "source_work_j", "ohmic_heat_j", "numerical_interface_heat_j",
                                                         "minimum_high_order_weight", "fluid_energy_balance_residual_j",
                                                         "ct_dissipation_heat_j", "response_energy_balance_residual_j", "dual_energy_change_j",
                                                         "dual_energy_absolute_change_j", "internal_energy_j",
                                                         "gas_energy_j", "induced_field_energy_j",
                                                         "dual_energy_cold_cell_steps", "dual_energy_cold_discrepancy_j",
                                                         "dual_energy_disagreement_cell_steps", "step_retries",
                                                         "floor_mass_kg", "floor_energy_j", "floor_cell_steps", "mass_kg",
                                                         "outflow_mass_kg", "outward_fluid_energy_j", "open_line_loss_j",
                                                         "response_energy_balance_relative")})
        history.append(row)
        z, rs, psi = separatrix(model)
        q = model._primitive(model.state, model.time)
        snapshots[f"t{len(history)-1:02d}"] = np.concatenate((psi.ravel(), q[..., 0].ravel(), q[..., 3].ravel()))
        velocity[f"t{len(history)-1:02d}"] = np.concatenate((q[..., 1].ravel(), q[..., 2].ravel()))
        print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
        (out/"history.json").write_text(json.dumps(history, indent=1))
        np.savez_compressed(out/"snapshots.npz", rf=model.mesh.rf, zf=model.mesh.zf, **snapshots)
        np.savez_compressed(out/"velocity.npz", **velocity)
    return history


def parse(field, text):
    kind = type(getattr(Nucte(), field))
    if kind is bool:
        if text not in ("True", "False"):
            raise ValueError(f"{field} must be True or False")
        return text == "True"
    return kind(text)


if __name__ == "__main__":
    args = dict(a.split("=") for a in sys.argv[2:])
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path("outputs/nucte_formation"),
         **{k: parse(k, v) for k, v in args.items()})
