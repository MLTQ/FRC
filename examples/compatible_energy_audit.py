"""Independently decompose accepted shock steps into spatial and RK energy defects."""
from pathlib import Path
import argparse
import hashlib
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.conservative_benchmarks import from_split_initial
from plasma_bench.shock_benchmarks import standing_shock
from plasma_bench.magnetic_quadrature import stiffness_action
from .conservative_measurements import provenance
from .verification_records import write_record


def quadratic_step_defect(solver, stages, dt):
    """Quadratic invariant defect from independent RK coefficient algebra.

    Constant external bias only. Recomplete each flux-rate direction through
    the vacuum and contract its full-space minus fluid-region magnetic form.
    The actual energy change is not used to construct this prediction.
    """
    m, mesh = solver.mapping, solver.mesh
    count = len(stages)
    if count == 2:
        weights = np.array([.5, .5]); matrix = np.array([[0., 0.], [1., 0.]])
    elif count == 3:
        weights = np.array([1/6, 1/6, 2/3])
        matrix = np.array([[0., 0., 0.], [1., 0., 0.], [.25, .25, 0.]])
    else:
        raise ValueError("Expected RK2 or RK3 stage rates")
    full_rates, fluid_actions, defects = [], [], []
    for state, time, total, rhs in stages:
        assert not np.any(solver.sources.at(time).dpsi_dt)
        du, dpsi, pout = rhs[:3]
        data = np.zeros(m.field.shape)
        data[:mesh.nr+1, m.z_offset:m.z_offset+mesh.nz+1] = dpsi
        full_rates.append(m.interior(m.harmonic_extension(data)))
        fluid_actions.append(stiffness_action(mesh, dpsi))
        ku = m.k@m.interior(state.induced_psi)
        weak_load = m.fluid_flux(m.full(ku))-stiffness_action(mesh,total)
        defects.append(-pout+2*np.pi/MU0*np.sum(weak_load*dpsi))
    gram = np.empty((count,count))
    for i in range(count):
        for j in range(count):
            gram[i,j] = np.pi/MU0*(full_rates[i]@(m.k@full_rates[j])-
                                  np.sum(stages[i][3][1]*fluid_actions[j]))
    coefficients = np.outer(weights,weights)-weights[:,None]*matrix-weights[None,:]*matrix.T
    spatial = dt*float(weights@defects)
    temporal = dt*dt*float(np.sum(coefficients*gram))
    return {"spatial_prediction_j": spatial, "temporal_prediction_j": temporal,
            "maximum_stage_spatial_defect_w": float(np.max(abs(np.asarray(defects)))),
            "gram_symmetry_error": float(np.max(abs(gram-gram.T)))}


def run(cells=96, order=3, cfl=.2, directory="outputs/compatible_energy_verification/time_audit"):
    solver = from_split_initial(standing_shock(cells,cfl=cfl)[0], energy_quadrature="staggered",
                                boundary_poynting="variational", time_order=order)
    original = solver._forward_euler
    stages = []
    def observe(state,time,total,dt):
        rhs, report = original(state,time,total,dt)
        stages.append((state,time,total,rhs))
        return rhs, report
    solver._forward_euler = observe
    directory = Path(directory); directory.mkdir(parents=True,exist_ok=True)
    path = directory/f"shock_{cells}_rk{order}_cfl{cfl}.json"
    record = {"radial_cells": cells, "time_order": order, "cfl": cfl,
              "provenance": provenance(), "audit_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "completed": False, "snapshots": []}
    accumulated_spatial = accumulated_temporal = max_stage = max_identity = 0.
    start_response = solver.initial["response_energy_j"]
    previous_response = start_response
    for end in (.005,.01,.015,.02):
        while solver.time < end:
            remaining = end-solver.time
            if remaining <= 4*np.finfo(float).eps*max(abs(end),1e-300):
                solver.time = end; break
            stages.clear()
            dt = solver.step(min(solver.timestep_limit(),remaining))
            audit = quadratic_step_defect(solver,stages,dt)
            response = solver._integrals(solver.state,solver.time)["response_energy_j"]
            measured = response-previous_response
            identity_error = measured-audit["spatial_prediction_j"]-audit["temporal_prediction_j"]
            max_identity = max(max_identity,abs(identity_error))
            max_stage = max(max_stage,audit["maximum_stage_spatial_defect_w"])
            accumulated_spatial += audit["spatial_prediction_j"]
            accumulated_temporal += audit["temporal_prediction_j"]
            previous_response = response
        record["snapshots"].append({"time_s": solver.time, "steps": solver.steps,
                                    "measured_response_change_j": response-start_response,
                                    "accumulated_spatial_prediction_j": accumulated_spatial,
                                    "accumulated_temporal_prediction_j": accumulated_temporal,
                                    "maximum_stage_spatial_defect_w": max_stage,
                                    "maximum_step_identity_error_j": max_identity})
        write_record(path,record)
        print(record["snapshots"][-1],flush=True)
    record.update(completed=True,diagnostics=solver.diagnostics())
    write_record(path,record)
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells",type=int,default=96)
    parser.add_argument("--order",type=int,default=3)
    parser.add_argument("--cfl",type=float,default=.2)
    args = parser.parse_args()
    run(args.cells,args.order,args.cfl)
