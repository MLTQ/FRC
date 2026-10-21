"""End-to-end checks of the rod-driven toroidal field in the 12" annular chamber (docs/EXPERIMENT_PLAYGROUND.md,
"Toroidal field"; companion rod_verification.md).

python -m examples.rod_verification OUT static [nr ...]     uniform B_z + B_theta = mu0 I/(2 pi r), uniform p: at rest?
python -m examples.rod_verification OUT vacuum [nr ...]     rod bank fired into a resistive gas: B_theta -> mu0 I(t)/(2 pi r)
                                                            with I(t) from the coupled circuit, against the plasma-free circuit
python -m examples.rod_verification OUT closure RUN ...     rod energy closure of saved formation runs (bank -> line, gas)

The gas, numerics and inner tube are those of frc12_inner_tube (formation parameters from nucte_formation.build).
Writes OUT/<case>.json and prints a table."""
import json
import sys
import time
from dataclasses import replace
from pathlib import Path
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.electromagnetics import VacuumField
from plasma_bench.mhd_state import from_primitive, with_toroidal
from plasma_bench.conservative_coupled import ConservativeVacuumMHD, ConservativeState
from plasma_bench.coil_circuit import simulate
from examples.nucte_formation import Nucte, build, inner_rows, parse
from examples.experiments import PRESETS, rod_bank

STATIC_CURRENT, STATIC_BIAS, STATIC_END = 40e3, -.05, 5e-6   # A, T, s (about 14 radial Alfven crossings at the wall)
VACUUM_ETA, VACUUM_FIRE, VACUUM_END = 5e-3, .5e-6, 20e-6      # ohm m (gas diffusion e-fold ~0.33 us), s, s


def case(nr, **extra):
    """frc12_inner_tube at nr x 6 nr with overrides."""
    return replace(Nucte(), **{**PRESETS["frc12_inner_tube"], "nr": nr, "nz": 6*nr, **extra})


def radius_squared(mesh):
    return .5*(mesh.rf[:-1]**2+mesh.rf[1:]**2)[:, None]


def static(nr, order=2):
    """Screw pinch at rest: uniform B_z (bias), the rod's vacuum B_theta = K/r (discrete g = K/<r^2>, zero discrete
    current) and uniform density and pressure (the 2 eV fill). Both fields are current-free, so the exact solution is
    rest for all time. Formation parameters; the threading r B_theta is held at K."""
    k = MU0*STATIC_CURRENT/(2*np.pi)
    c = case(nr, drive="ideal", pi_half_length=0., conduction_limit=0., peak_t=STATIC_BIAS, bias_t=STATIC_BIAS,
             rod_bank=rod_bank(40, -20e-6))
    formation, setup = build(c)                    # for its mesh, grid and parameters
    mesh, grid = formation.mesh, formation.mapping.field
    i0, _ = inner_rows(c)
    sources = VacuumField(grid, (), bias=STATIC_BIAS)
    psi = sources.at(0.).psi[i0:i0+c.nr+1, c.end_cells:c.end_cells+c.nz+1]
    state = from_primitive(mesh, psi, setup["density_kg_m3"], setup["pressure_pa"], energy_quadrature=c.energy_quadrature)
    r2 = radius_squared(mesh)
    conserved = with_toroidal(mesh, state.conserved, 0., k/r2)
    model = ConservativeVacuumMHD(mesh, sources, ConservativeState(conserved, np.zeros(grid.shape)),
                                  replace(formation.parameters, spatial_order=order, toroidal_boundary=lambda t: k))
    rho = setup["density_kg_m3"]
    b_in = np.hypot(STATIC_BIAS, k/mesh.rf[0])
    v_a = b_in/np.sqrt(MU0*rho)
    rows, start = [], time.time()
    for t in np.linspace(STATIC_END/5, STATIC_END, 5):
        model.advance(t)
        q = model._primitive(model.state, model.time)
        d = model.diagnostics()
        a = r2*q[..., 7]
        rows.append(dict(time_s=model.time, steps=model.steps, max_speed_m_s=float(np.hypot(q[..., 1], q[..., 2]).max()),
                         max_vtheta_m_s=float(abs(np.sqrt(r2)*q[..., 6]).max()),
                         max_rbtheta_deviation=float(abs(a/k-1).max()),
                         max_pressure_deviation=float(abs(q[..., 3]/setup["pressure_pa"]-1).max()),
                         kinetic_energy_j=float(np.sum(.5*q[..., 0]*(q[..., 1]**2+q[..., 2]**2)*mesh.volumes))
                         + d["toroidal_kinetic_energy_j"],
                         boundary_work_j=d["toroidal_boundary_work_j"],
                         response_energy_residual_j=d["response_energy_balance_residual_j"],
                         fluid_energy_residual_j=d["fluid_energy_balance_residual_j"]))
    return dict(case="static", nr=nr, nz=c.nz, spatial_order=order, rod_current_a=STATIC_CURRENT, bias_t=STATIC_BIAS,
                alfven_speed_inner_wall_m_s=float(v_a), toroidal_energy_j=model.initial["toroidal_magnetic_energy_j"],
                wall_s=time.time()-start, at=rows)


def rod_energy(row, rod, bank):
    """Rod energy budget [J] at a saved state: bank energy against capacitor, stray, crowbar, dissipation, the line's
    field outside the gas 1/2 (L - L_gas) I^2, the gas B_theta energy at t = 0 and the boundary Poynting work.
    L_gas is the discrete vacuum's (exactly the energy of g = K/<r^2>). I is the circuit branch current."""
    s = row["rod_circuit"]
    l_line, l_gas = rod["line_inductance_h"], rod["gas_vacuum_inductance_h"]
    i_m, i_c, v_c = s["module_current_a"][0], s["crowbar_current_a"][0], s["capacitor_voltage_v"][0]
    e0 = .5*bank.capacitance*bank.voltage**2
    terms = dict(capacitor_j=.5*bank.capacitance*v_c**2, stray_j=.5*bank.stray*i_m**2,
                 crowbar_j=.5*bank.crowbar_inductance*i_c**2, dissipated_j=s["dissipated_j"],
                 line_field_j=.5*(l_line-l_gas)*s["branch_current_a"]**2,
                 gas_initial_j=row["_gas_initial_j"], boundary_work_j=row.get("toroidal_boundary_work_j", 0.))
    return e0, terms, e0-sum(terms.values())


def vacuum(nr):
    """20 kA rod bank fired at 0.5 us into a uniformly resistive gas (eta 5e-3 ohm m, diffusion e-fold ~0.33 us) under a
    constant bias field: the gas cannot hold toroidal flux, so r B_theta -> mu0 I(t)/(2 pi) everywhere and I(t) follows
    the plasma-free circuit (coil_circuit.simulate, an independent adaptive integrator)."""
    c = case(nr, drive="ideal", pi_half_length=0., conduction_limit=0., peak_t=-.05, bias_t=-.05, crowbar_decay_s=1.,
             resistivity_model="uniform", resistivity=VACUUM_ETA, chodura=0., end_s=VACUUM_END,
             rod_bank=rod_bank(20, VACUUM_FIRE))
    model, setup = build(c)
    rod, bank = model.rod, c.rod_bank
    ref = simulate([[rod.inductance]], rod.modules, VACUUM_END, samples=200001)
    mesh = model.mesh
    r2 = radius_squared(mesh)
    mid = c.nz//2
    rows, start = [], time.time()
    for t in np.arange(2e-6, VACUUM_END+1e-12, 2e-6):
        model.advance(t, on_step=model.coupler)
        a = r2*model.state.conserved[..., 5]
        k = model._linked(model.time)
        i_ref = float(np.interp(model.time, ref.time, ref.section_current[:, 0]))
        d = model.diagnostics()
        row = dict(time_s=model.time, steps=model.steps, rod_circuit=model.rod.snapshot(model),
                   toroidal_boundary_work_j=d["toroidal_boundary_work_j"], _gas_initial_j=0.)
        e0, terms, residual = rod_energy(row, setup["rod"], bank)
        # B_theta at r = 8 cm on the midplane from the cell r B_theta (B_theta = a/r)
        bt8 = float(np.interp(.08, np.sqrt(r2[:, 0]), .5*(a[:, mid-1]+a[:, mid]))/.08)
        rows.append(dict(time_s=model.time, steps=model.steps, rod_current_a=row["rod_circuit"]["rod_current_a"],
                         plasma_free_current_a=i_ref,
                         current_relative_difference=(row["rod_circuit"]["rod_current_a"]-i_ref)/max(abs(i_ref), 1.),
                         max_rbtheta_over_k_minus_1=float(abs(a/k-1).max()) if k else None,
                         midplane_rbtheta_over_k_minus_1=float(abs(a[:, [mid-1, mid]]/k-1).max()) if k else None,
                         btheta_8cm_t=bt8, vacuum_btheta_8cm_t=MU0*i_ref/(2*np.pi*.08),
                         plasma_linkage_wb=row["rod_circuit"]["plasma_linkage_wb"],
                         boundary_work_j=terms["boundary_work_j"], toroidal_magnetic_energy_j=d["toroidal_magnetic_energy_j"],
                         ohmic_heat_j=d["ohmic_heat_j"], rod_energy_residual_j=residual, rod_bank_energy_j=e0,
                         response_energy_residual_j=d["response_energy_balance_residual_j"],
                         max_speed_m_s=d["maximum_speed_m_s"]))
    return dict(case="vacuum", nr=nr, nz=c.nz, resistivity_ohm_m=VACUUM_ETA, rod=setup["rod"], wall_s=time.time()-start,
                at=rows)


def closure(run):
    """Rod energy closure of a saved formation run, per saved state (history rows)."""
    run = Path(run)
    setup = json.loads((run/"setup.json").read_text())
    bank = parse("rod_bank", json.dumps(setup["case"]["rod_bank"]))
    rod = setup["rod"]
    history = json.loads((run/"history.json").read_text())
    i0 = history[0]["rod_circuit"]["rod_current_a"]
    gas0 = .5*rod["gas_vacuum_inductance_h"]*i0**2
    rows, worst = [], 0.
    for h in history[1:]:
        e0, terms, residual = rod_energy({**h, "_gas_initial_j": gas0}, rod, bank)
        worst = max(worst, abs(residual))
        rows.append(dict(time_s=h["time_s"], residual_j=residual, residual_over_bank=residual/e0,
                         residual_over_gas_energy=residual/max(terms["boundary_work_j"]+gas0, 1e-30), **terms))
    return dict(case="closure", run=str(run), bank_energy_j=.5*bank.capacitance*bank.voltage**2, gas_initial_j=gas0,
                worst_residual_j=worst, end_s=history[-1]["time_s"], at=rows[::max(1, len(rows)//40)]+rows[-1:])


def main(out, kind, args):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    if kind == "closure":
        results = [closure(run) for run in args]
        for r in results:
            last = r["at"][-1]
            print(f"{r['run']}: bank {r['bank_energy_j']:.1f} J, to {r['end_s']*1e6:.1f} us; worst |residual| "
                  f"{r['worst_residual_j']:.3g} J ({r['worst_residual_j']/r['bank_energy_j']:.2e} of the bank); at the end: "
                  + ", ".join(f"{k} {v:.4g}" for k, v in last.items() if k.endswith("_j")))
    else:
        results = [dict(static=static, static_first=lambda n: static(n, order=1), vacuum=vacuum)[kind](int(n)) for n in (args or ["12", "24"])]
        for r in results:
            print(json.dumps({k: v for k, v in r.items() if k != "at"}))
            for row in r["at"]:
                print("  "+json.dumps({k: (float(f"{v:.6g}") if isinstance(v, float) else v) for k, v in row.items()
                                        if not isinstance(v, dict)}))
    (out/f"{kind}.json").write_text(json.dumps(results, indent=1))
    return results


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
