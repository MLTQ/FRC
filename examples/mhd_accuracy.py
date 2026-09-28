"""Independent time-refinement and first-order baseline comparisons."""
from dataclasses import asdict
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.mhd import AxisymmetricMHD, MHDParameters
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.mhd_state import primitive
from plasma_bench.mhd_benchmarks import radial_fast_mode, advected_flux, resistive_bessel_mode, sod_initial, sod_exact


def relative_l2(actual, expected, weight=1):
    return float(np.sqrt(np.sum((actual-expected)**2*weight)/np.sum(expected**2*weight)))


def temporal_refinement(*, smooth):
    mesh = MHDMesh(1, 1, 24, 4)
    epsilon = 1e-3 if smooth else 1e-5
    initial, omega = radial_fast_mode(mesh, epsilon=epsilon)
    end_time = .04 if smooth else np.pi/(2*omega)
    baseline = AxisymmetricMHD(mesh, initial).timestep_limit()
    counts = (16, 32, 64, 512) if smooth else [int(np.ceil(end_time/(f*baseline))) for f in (.8, .4, .2, .025)]
    velocities = []
    for count in counts:
        solver = AxisymmetricMHD(mesh, initial)
        for _ in range(count):
            solver.step(end_time/count)
        velocities.append(primitive(mesh, solver.state, 5/3)[..., 1])
    return [{"steps": count, "reference_steps": counts[-1], "end_time_s": end_time,
             "epsilon": epsilon, "mesh": asdict(mesh), "parameters": asdict(solver.parameters),
             "velocity_relative_L2_vs_time_reference": relative_l2(velocity, velocities[-1])}
            for count, velocity in zip(counts[:-1], velocities[:-1])]


def reference_comparison(metrics):
    cases = []
    for name, count, field in (("fast_wave", 96, "velocity_relative_L2"),
                               ("advection", 48, "flux_relative_L2"),
                               ("resistive_decay", 64, "flux_relative_L2"),
                               ("sod", 160, "density_L1")):
        key = "axial_cells" if name == "sod" else "radial_cells"
        row = next(v for v in metrics[name] if v[key] == count)
        mesh = MHDMesh(**row["mesh"])
        par = MHDParameters(**{**row["parameters"], "spatial_order": 1})
        if name == "fast_wave":
            initial, omega = radial_fast_mode(mesh)
            expected, _ = radial_fast_mode(mesh, time=np.pi/(2*omega))
        elif name == "advection":
            initial, expected = advected_flux(mesh), advected_flux(mesh, time=.15)
        elif name == "resistive_decay":
            initial, k = resistive_bessel_mode(mesh)
            expected = initial.psi*np.exp(-.5*k*k*.01)
        else:
            initial = sod_initial(mesh)
            expected = sod_exact(mesh.z, .12)[0]
        solver = AxisymmetricMHD(mesh, initial, par)
        solver.advance(row["time_s"])
        if name == "fast_wave":
            error = relative_l2(primitive(mesh, solver.state, par.gamma)[..., 1], primitive(mesh, expected, par.gamma)[..., 1], mesh.volumes)
        elif name == "advection":
            error = relative_l2(solver.state.psi, expected.psi)
        elif name == "resistive_decay":
            error = relative_l2(solver.state.psi, expected)
        else:
            select = abs(mesh.z) < .4
            error = float(np.mean(abs(primitive(mesh, solver.state, par.gamma)[1, select, 0]-expected[select])))
        cases.append({"case": name, "cells": count, "quantity": field, "first_order_error": error,
                      "second_order_error": row[field], "error_reduction_factor": error/row[field],
                      "mesh": asdict(mesh), "reference_parameters": asdict(par),
                      "reference_diagnostics": solver.diagnostics()})
    return cases
