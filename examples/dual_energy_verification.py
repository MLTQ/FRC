"""Dual-energy and trace-bound checks on the existing conservative fixtures: warm flow, shock accuracy, low-beta adiabaticity."""
import json
import sys
import time
from pathlib import Path
import numpy as np
from plasma_bench.conservative_benchmarks import from_split_initial
from plasma_bench.coupled_benchmarks import coil_pulse
from plasma_bench.coupled_hydro import gas_from_primitive
from plasma_bench.shock_benchmarks import standing_shock
from plasma_bench.mhd_state import NonphysicalState
from examples.conservative_measurements import shock_errors

KEEP = ("time_s", "steps", "response_energy_balance_relative", "fluid_energy_balance_residual_j", "source_work_j",
        "fluid_total_energy_j", "limited_euler_stages", "dual_energy_change_j", "dual_energy_cold_discrepancy_j",
        "dual_energy_cold_cell_steps", "dual_energy_disagreement_cell_steps", "dual_energy_max_cold_ratio")
# (label, overrides); the default dual-energy configuration is pressure mode, threshold 0.1, detector 0.2.
CONFIGS = (("off", dict(dual_energy=False)), ("default", dict(dual_energy=True)),
           ("no_detector", dict(dual_energy=True, dual_energy_shock_threshold=None)),
           ("detector_0.05", dict(dual_energy=True, dual_energy_shock_threshold=.05)),
           ("reset_mode", dict(dual_energy=True, dual_energy_mode="reset")),
           ("threshold_1e-3", dict(dual_energy=True, dual_energy_threshold=1e-3)),
           ("bounded_traces", dict(local_trace_bounds=True)),
           ("dual_bounded", dict(dual_energy=True, local_trace_bounds=True)),
           ("entropy_reset", dict(dual_energy=True, dual_energy_mode="reset", dual_energy_criterion="entropy")),
           ("ct_flux", dict(dual_energy=True, dual_energy_mode="reset", dual_energy_criterion="entropy",
                            magnetic_energy_flux="ct")),
           ("ct_compatible", dict(dual_energy=True, dual_energy_mode="reset", dual_energy_criterion="entropy",
                                  magnetic_energy_flux="ct", energy_quadrature="staggered",
                                  boundary_poynting="variational", time_order=3)),
           ("compatible", dict(dual_energy=True, dual_energy_mode="reset", dual_energy_criterion="entropy",
                               energy_quadrature="staggered", boundary_poynting="variational", time_order=3)))


def run(solver, end, every):
    start, failure, invariant = time.time(), None, []
    try:
        while solver.time < end:
            solver.advance(min(end, solver.time+every))
            yield_state = solver._primitive(solver.state, solver.time)
            invariant.append(yield_state)
    except (NonphysicalState, ValueError, RuntimeError) as error:
        failure = f"{type(error).__name__}: {error}"
    d = solver.diagnostics()
    return dict(completed=failure is None, failure=failure, wall_s=time.time()-start, **{k: d[k] for k in KEEP}), invariant


def shock(overrides, cells):
    split, reference = standing_shock(cells)
    solver = from_split_initial(split, **overrides)
    row, _ = run(solver, .02, .02)
    if row["completed"]:
        row["shock_errors"] = shock_errors(solver, reference)[0]
    return row, solver.state


def pulse(overrides, cells):
    solver = from_split_initial(coil_pulse(cells, vacuum_boundary="multipole"), **overrides)
    row, _ = run(solver, .1, .1)
    return row, solver.state


def low_beta(overrides, cells, pressure=1e-6):
    """Coil pulse into beta ~ 2e-6 gas, no physical or interface resistive heating.

    Independent reference: smooth adiabatic flow keeps A = p/rho^gamma at its uniform initial value.
    """
    split = coil_pulse(cells, vacuum_boundary="multipole")
    split.state.gas = gas_from_primitive(split.mesh, 1., pressure)
    solver = from_split_initial(split, resistivity=0., interface_diffusion=False, **overrides)
    row, states = run(solver, .1, .01)
    gamma = solver.parameters.gamma
    row["max_adiabatic_invariant_error"] = max((float(abs(q[..., 3]/q[..., 0]**gamma/pressure-1).max()) for q in states),
                                               default=None)
    return row


def main(out=Path("outputs/dual_energy_verification"), cells=16):
    out.mkdir(parents=True, exist_ok=True)
    result = {}
    for label, overrides in CONFIGS:
        s_row, s_state = shock(overrides, cells)
        p_row, p_state = pulse(overrides, cells)
        result[label] = dict(shock=s_row, pulse=p_row, low_beta=low_beta(overrides, cells), states=(s_state, p_state))
    off_states = result["off"].pop("states")
    for label in result:
        states = result[label].pop("states", off_states)
        result[label]["pulse_max_conserved_difference_from_off"] = float(abs(states[1].conserved-off_states[1].conserved).max())
    (out/"results.json").write_text(json.dumps(result, indent=1))
    for label, r in result.items():
        e = r["shock"].get("shock_errors", {})
        lb = r["low_beta"]
        print(f"{label:15s} shock {'ok' if r['shock']['completed'] else 'FAIL'} rhoL1={e.get('density_L1_kg_m3', float('nan')):.5f} "
              f"pL1={e.get('pressure_L1_pa', float('nan')):.4f} | pulse diff from off={r['pulse_max_conserved_difference_from_off']:.3g} "
              f"| low-beta {'ok' if lb['completed'] else 'FAIL t=%.4f' % lb['time_s']} invariant={lb['max_adiabatic_invariant_error'] or float('nan'):.3g} "
              f"fluid residual={lb['fluid_energy_balance_residual_j']:.2g} J reset dE={lb['dual_energy_change_j']:.2g} J "
              f"cold discrepancy={lb['dual_energy_cold_discrepancy_j']:.2g} J")


if __name__ == "__main__" and "--floor" not in sys.argv:
    main(*(sys.argv[1:2] and [Path(sys.argv[1])]))


def floor_accounting(cells=16):
    """Declared minimum-density closure: inactive -> bit-identical; active -> balances close with its counted terms."""
    out = {}
    base = from_split_initial(coil_pulse(cells, vacuum_boundary="multipole"), dual_energy=True)
    base.advance(.05)
    for label, floor in (("inactive", 1e-6), ("active", 1.05)):
        s = from_split_initial(coil_pulse(cells, vacuum_boundary="multipole"), dual_energy=True, minimum_density=floor,
                               floor_temperature_ev=1e-3, floor_ion_mass=1.)
        s.advance(.05)
        d = s.diagnostics()
        q = s._primitive(s.state, s.time)
        out[label] = dict(max_abs_difference_from_no_floor=float(abs(s.state.conserved-base.state.conserved).max()),
                          minimum_density=float(q[..., 0].min()), floor_mass_kg=d["floor_mass_kg"],
                          floor_energy_j=d["floor_energy_j"], floor_cell_steps=d["floor_cell_steps"],
                          mass_relative_change_excluding_floor=d["mass_relative_change"],
                          fluid_energy_residual_j=d["fluid_energy_balance_residual_j"])
    return out


if __name__ == "__main__" and "--floor" in sys.argv:
    print(json.dumps(floor_accounting(), indent=1))
