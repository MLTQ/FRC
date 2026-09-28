"""Reproducible numerical measurements for the fixed-chamber coupling reference."""
from dataclasses import asdict
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.coupled_benchmarks import (aligned_meshes, compact_flux, compact_force,
    cell_average, passive_advection, passive_diffusion, gaussian_diffusion, coil_pulse)
from plasma_bench.coupling_mesh import CouplingMesh
from plasma_bench.coupled_magnetic import MagneticStepper
from plasma_bench.coupled_hydro import gas_from_primitive
from plasma_bench.coupled import CoupledState, VacuumCoupledMHD
from plasma_bench.electromagnetics import VacuumField
from plasma_bench.vacuum_benchmarks import PROBES, infinite_flux, finite_box_flux
from plasma_bench.waveforms import to_dict


def configuration(solver):
    return {"fluid_mesh": asdict(solver.mesh), "field_grid": asdict(solver.mapping.field),
            "parameters": asdict(solver.parameters), "bias_t": solver.sources.bias,
            "coils": [{**asdict(coil), "waveform": to_dict(coil.waveform)} for coil in solver.sources.coils]}


def force_measurements():
    rows = []
    for cells in (8, 16, 32, 64):
        mesh, grid = aligned_meshes(cells)
        mapping = CouplingMesh(mesh, grid)
        u = compact_flux(grid.r[:, None], grid.z[None, :])
        total = u+.5*np.sqrt(MU0)*grid.r[:, None]**2
        force = MagneticStepper(mapping, 0).advection(total).T@((2*np.pi/MU0)*(mapping.k@mapping.interior(u)))
        measured = [f.reshape(mesh.nr, mesh.nz)/mesh.volumes for f in np.split(force, 2)]
        exact = cell_average(mesh, compact_force)
        err = sum((a-b)**2 for a, b in zip(measured, exact))
        norm = sum(a*a for a in exact)
        rows.append({"radial_cells": cells, "relative_volume_L2": float(np.sqrt(np.sum(err*mesh.volumes)/np.sum(norm*mesh.volumes))),
                     "maximum_axis_force_error_n_m3": float(np.sqrt(err[0]).max()),
                     "fluid_mesh": asdict(mesh), "field_grid": asdict(grid),
                     "amplitude": .01, "width_m": .65, "bias_t": float(np.sqrt(MU0))})
    return rows


def passive_measurements(kind):
    rows = []
    for cells in (12, 24, 48):
        solver = passive_advection(cells) if kind == "advection" else passive_diffusion(cells)
        end = .05 if kind == "advection" else .03
        solver.advance(end)
        grid = solver.mapping.field
        r, z = grid.r[:, None], grid.z[None, :]
        exact = compact_flux(r, z, width=.5, amplitude=1e-6, shift=.025) if kind == "advection" else gaussian_diffusion(r, z, end)
        width = .8 if kind == "advection" else .7
        select = (r < width)&(abs(z) < width)
        error = float(np.linalg.norm((solver.state.induced_psi-exact)[select])/np.linalg.norm(exact[select]))
        rows.append({"radial_cells": cells, "relative_flux_L2": error, "comparison_window_m": width,
                     "configuration": configuration(solver), "diagnostics": solver.diagnostics()})
    return rows


def temporal_measurements():
    mesh, grid = aligned_meshes(8)
    u = compact_flux(grid.r[:, None], grid.z[None, :], amplitude=.3)
    results, diagnostics = [], []
    for steps in (8, 16, 32, 256):
        solver = VacuumCoupledMHD(mesh, VacuumField(grid, (), bias=np.sqrt(MU0)),
                                 CoupledState(gas_from_primitive(mesh, 1., 1.), u))
        for _ in range(steps):
            solver.step(.02/steps)
        results.append(solver.state.gas[..., 1:3])
        diagnostics.append(solver.diagnostics())
    errors = [float(np.linalg.norm(v-results[-1])) for v in results[:-1]]
    return {"step_counts": [8, 16, 32], "reference_step_count": 256, "end_time_s": .02,
            "momentum_errors": errors, "ratios": [errors[0]/errors[1], errors[1]/errors[2]],
            "initial_flux_amplitude": .3, "initial_flux_width_m": .65,
            "configuration": configuration(solver), "diagnostics": diagnostics}


def vacuum_measurements():
    exact = infinite_flux()
    rows = []
    for cells, extent, exact_boundary in ((8, 2, True), (16, 2, True), (32, 2, True),
                                         (32, 2, False), (32, 4, False), (32, 8, False)):
        value, residual = finite_box_flux(cells, extent, exact_boundary=exact_boundary)
        rows.append({"cells_per_m": cells, "outer_radius_and_half_length_m": extent,
                     "exact_boundary": exact_boundary, "relative_probe_L2": float(np.linalg.norm(value-exact)/np.linalg.norm(exact)),
                     "relative_linear_residual": residual, "flux_wb_per_rad": value.tolist()})
    return {"probes_r_z_m": PROBES.tolist(), "infinite_flux_wb_per_rad": exact.tolist(),
            "quadrature_order": 32, "quadrature_32_vs_48_relative": float(np.linalg.norm(infinite_flux(order=48)-exact)/np.linalg.norm(exact)),
            "cases": rows}


def pulse_measurement(cells, extent=2, cfl=.2, *, snapshots=False):
    solver = coil_pulse(cells, extent=extent, cfl=cfl)
    history = [solver.diagnostics()]
    for t in np.linspace(0, .1, 21)[1:]:
        solver.advance(float(t))
        history.append(solver.diagnostics())
    row = {"radial_cells": cells, "extent": extent, "configuration": configuration(solver),
           "diagnostics": solver.diagnostics(), "history": history}
    arrays = {}
    if snapshots:
        current = np.divide(solver.mapping.k@solver.mapping.interior(solver.state.induced_psi),
                            MU0*solver.mapping.node_area, out=np.full(solver.mapping.node_count, np.nan),
                            where=solver.mapping.active)
        arrays = {"gas": solver.state.gas, "induced_psi": solver.state.induced_psi,
                  "total_psi": solver.total_flux(), "current": current.reshape(solver.mapping.field.nr-2, solver.mapping.field.nz-2),
                  "fluid_r_faces": solver.mesh.rf, "fluid_z_faces": solver.mesh.zf,
                  "field_r": solver.mapping.field.r, "field_z": solver.mapping.field.z}
    return row, arrays


def all_measurements():
    metrics = {"force": force_measurements()}
    for kind in ("advection", "diffusion"):
        metrics[kind] = passive_measurements(kind)
        print(f"Completed {kind} refinement", flush=True)
    metrics["temporal"] = temporal_measurements()
    metrics["vacuum"] = vacuum_measurements()
    metrics["pulse"] = []
    arrays = {}
    for cells, extent, cfl in ((8, 2, .2), (16, 2, .2), (24, 2, .2), (32, 2, .2),
                               (16, 3, .2), (16, 4, .2), (16, 2, .1)):
        row, fields = pulse_measurement(cells, extent, cfl, snapshots=cells == 32)
        metrics["pulse"].append(row)
        arrays.update(fields)
        print(f"Completed pulse: {cells} radial cells, extent={extent}, CFL={cfl}", flush=True)
    return metrics, arrays
