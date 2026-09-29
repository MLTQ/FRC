"""Differentiate response energy independently to audit shared interface rates."""
from pathlib import Path
from dataclasses import asdict
import argparse
import hashlib
import json
import numpy as np
from scipy.constants import mu_0 as MU0
from .conservative_measurements import provenance
from plasma_bench.shock_benchmarks import standing_shock
from plasma_bench.conservative_benchmarks import from_split_initial


def rate_parts(solver, rhs):
    mesh, mapping = solver.mesh, solver.mapping
    br, bz = mesh.cell_fields(mapping.fluid_flux(solver.total_flux()))
    dbr, dbz = mesh.cell_fields(rhs[1])
    cell_magnetic = float(np.sum((br*dbr+bz*dbz)*mesh.volumes)/MU0)
    field_rate = np.zeros(mapping.field.shape)
    field_rate[:mesh.nr+1, mapping.z_offset:mapping.z_offset+mesh.nz+1] = rhs[1]
    ku = mapping.k@mapping.interior(solver.state.induced_psi)
    nodal = float(2*np.pi/MU0*(ku@mapping.grid_vector(field_rate)))
    return {"outward_fluid_power_w": rhs[2], "cell_magnetic_power_w": cell_magnetic,
            "nodal_total_flux_work_w": nodal, "response_defect_rate_w": -rhs[2]-cell_magnetic+nodal}


def snapshot(solver):
    state, time = solver.state, solver.time
    total = solver.mapping.fluid_flux(solver.total_flux())
    high, _ = solver._fluxes(state, time)
    low, _ = solver._fluxes(state, time, spatial_order=1)
    dt = solver.timestep_limit()
    selected, report = solver._forward_euler(state, time, total, dt)
    rates = {name: rate_parts(solver, rhs) for name, rhs in
             (("high", high.rhs(solver.mesh)), ("low", low.rhs(solver.mesh)), ("selected", selected))}
    differences = []
    # The shock's prescribed field is constant. Completing each perturbed
    # vacuum differentiates the independently integrated full-space energy.
    assert not np.any(solver.sources.at(time).dpsi_dt)
    for fraction in (.03, .01, .003):
        epsilon = fraction*dt
        plus = solver._new_state(state.conserved+epsilon*selected[0], total+epsilon*selected[1], time)
        minus = solver._new_state(state.conserved-epsilon*selected[0], total-epsilon*selected[1], time)
        fd = (solver._integrals(plus, time)["response_energy_j"]-
              solver._integrals(minus, time)["response_energy_j"])/(2*epsilon)
        differences.append({"perturbation_s": epsilon, "derivative_w": fd,
                            "absolute_identity_error_w": abs(fd-rates["selected"]["response_defect_rate_w"])})
    return {"time_s": time, "timestep_limit_s": dt, "rates": rates,
            "selected_minus_high_w": rates["selected"]["response_defect_rate_w"]-rates["high"]["response_defect_rate_w"],
            "affected_cells": int(report.affected.sum()), "troubled_cells": int(report.troubled.sum()),
            "finite_differences": differences, "diagnostics": solver.diagnostics()}


def run(cells, directory="outputs/local_limiter_verification/energy_audit"):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    solver = from_split_initial(standing_shock(cells)[0])
    record = {"radial_cells": cells, "parameters": asdict(solver.parameters), "provenance": provenance(),
              "audit_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "snapshots": []}
    for time in np.linspace(0, .02, 5):
        solver.advance(float(time))
        row = snapshot(solver); record["snapshots"].append(row)
        # This audit adds no state changes and never resets a budget.
        print(cells, time, row["rates"], "selected-high", row["selected_minus_high_w"], flush=True)
        (directory/f"shock_{cells}.json").write_text(json.dumps(record, indent=2, allow_nan=False)+"\n")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, nargs="+", default=[96,192])
    args = parser.parse_args()
    for n in args.cells:
        run(n)
