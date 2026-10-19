"""Coarse NUCTE-III field-reversed theta-pinch formation on the fixed-chamber conservative MHD path."""
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

D_MASS = 2.01410*AMU


@dataclass(frozen=True)
class Nucte:
    """Published NUCTE-III geometry/fill (Takahashi 2000; Kuroyanagi 2014) plus explicit assumptions."""
    wall_radius: float = .115         # quartz tube, 23 cm ID, 2.0 m long
    half_length: float = 1.0
    coil_radius: float = .17          # 18 single-turn straps, 5 cm wide
    end_radius: float = .162          # 5 straps each end; (0.17/0.162)^2 = 1.10, case 1 mirror ratio
    strap_width: float = .05
    centre_straps: int = 18
    end_straps: int = 5
    filaments: int = 4                # per strap, equal share of the strap current
    fill_pa: float = 1.33             # D2
    gas_k: float = 300.
    pi_temperature_ev: float = 2.     # ASSUMED complete pre-ionization, Te = Ti
    pi_half_length: float = 0.        # >0: ionized only for |z| < this (theta-PI under the coil); beyond, floor-level plasma
    pi_edge_width: float = .03        # tanh transition width of the pre-ionized region [m]
    bias_t: float = -.06              # ASSUMED; published maximum is 0.09 T
    peak_t: float = .4                # ASSUMED vacuum centre field at end of rise
    quarter_period_s: float = 4e-6    # published rise time
    crowbar_decay_s: float = 100e-6   # ASSUMED L/R after crowbar
    resistivity: float = 1e-5         # ASSUMED uniform, Ohm m (Spitzer at ~13 eV); used when resistivity_model="uniform"
    resistivity_model: str = "uniform"  # or "spitzer": NRL eta_perp of local Te (Te=Ti, lnL=10, Te>=1 eV for eta)
    chodura: float = 0.               # Chodura anomalous-resistivity coefficient (standard ~0.1; 0 disables)
    nr: int = 16
    nz: int = 128
    exterior_cells: int = 6           # fixed induced-flux boundary just inside the coil
    end_cells: int = 2                # vacuum field cells beyond each tube end before the fixed boundary
    end_s: float = 15e-6
    cfl: float = .2
    dual_energy: bool = True          # entropy pressure in cold/low-beta cells (pressure mode, energy exact)
    max_step_halvings: int = 6        # timestep control: retry inadmissible steps with dt/2
    interface_diffusion: bool = True  # solver's declared numerical resistivity at the chamber boundary
    interface_speed_cap: float = 0.   # >0 caps the signal speed in that numerical resistivity [m/s]
    open_ends: bool = False           # outflow-only tube ends (open field-line end losses) instead of reflecting
    conduction_limit: float = 0.      # >0: flux-limited parallel conduction loss on open lines (f ~ 0.1)
    local_trace_bounds: bool = True   # reconstructed density/pressure bounded by neighbour minima
    floor_fraction: float = .01       # DECLARED gap closure: minimum density as a fraction of the fill (0 disables)


def strap_layout(c):
    w = c.strap_width
    centre = [(c.coil_radius, (k-(c.centre_straps-1)/2)*w) for k in range(c.centre_straps)]
    edge = c.centre_straps*w/2
    return centre+[(c.end_radius, s*(edge+(k+.5)*w)) for s in (-1, 1) for k in range(c.end_straps)]


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


def build(c):
    fluid = MHDMesh(c.wall_radius, c.half_length, c.nr, c.nz)
    grid = Grid(c.wall_radius+c.exterior_cells*fluid.dr, c.half_length+c.end_cells*fluid.dz,
                c.nr+c.exterior_cells+1, c.nz+2*c.end_cells+1)
    fil, current = equal_flux_currents(c)
    wave = centre_field(c)
    coils = tuple(Coil(r, z, wave, turns=float(i)/c.filaments, name=f"strap{k}")
                  for k, (strap, i) in enumerate(zip(fil, current)) for r, z in strap)
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
                                    interface_diffusion=c.interface_diffusion,
                                    interface_speed_cap=c.interface_speed_cap if c.interface_speed_cap > 0 else None,
                                    open_ends=c.open_ends, open_line_flux_limit=c.conduction_limit,
                                    open_line_wall_temperature_ev=c.pi_temperature_ev,
                                    local_trace_bounds=c.local_trace_bounds,
                                    minimum_density=c.floor_fraction*2*c.fill_pa/(KB*c.gas_k)*D_MASS if c.floor_fraction > 0 else None,
                                    floor_temperature_ev=c.pi_temperature_ev, floor_ion_mass=D_MASS)
    # Initial plasma is current-free in the bias: induced flux zero, total flux = source flux.
    probe = ConservativeState(np.zeros((c.nr, c.nz, 4)), np.zeros((grid.nr, grid.nz)))
    total = sources.at(0.).psi[:c.nr+1, c.end_cells:c.end_cells+c.nz+1]
    state = from_primitive(fluid, total, rho, p)
    model = ConservativeVacuumMHD(fluid, sources, ConservativeState(state.conserved, probe.induced_psi), params)
    return model, dict(fill_deuterons_m3=n_d, density_kg_m3=n_d*D_MASS, pressure_pa=2*n_d*c.pi_temperature_ev*QE,
                       strap_currents_a_per_t=current.tolist(),
                       centre_strap_share=float(current[:c.centre_straps].sum()/current.sum()))


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
                max_pressure_pa=float(q[..., 3].max()), max_speed_m_s=float(np.hypot(q[..., 1], q[..., 2]).max()))


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
            model.advance(target)
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
                                                         "dual_energy_cold_cell_steps", "dual_energy_cold_discrepancy_j",
                                                         "dual_energy_disagreement_cell_steps", "step_retries",
                                                         "floor_mass_kg", "floor_energy_j", "floor_cell_steps", "mass_kg",
                                                         "outflow_mass_kg", "outward_fluid_energy_j", "open_line_loss_j")})
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
