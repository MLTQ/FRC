"""Independent static and dynamic measurements of the infinite magnetic exterior."""
from dataclasses import asdict
import numpy as np
from plasma_bench.open_benchmarks import CurrentPacket, static_packet, nodal_probes, INTERIOR_PROBES, EXTERIOR_PROBES
from plasma_bench.open_boundary import multipoles
from plasma_bench.coupled_benchmarks import coil_pulse
from .coupled_measurements import configuration


def normal_derivative_jump(mapping, unknowns):
    grid, boundary = mapping.field, mapping.open_boundary
    full = mapping.full(unknowns)
    pr = np.gradient(full, grid.dr, axis=0, edge_order=2)[boundary.bi, boundary.bj]
    pz = np.gradient(full, grid.dz, axis=1, edge_order=2)[boundary.bi, boundary.bj]
    _, dr, dz = multipoles(grid.r[boundary.bi], grid.z[boundary.bj], boundary.scale, boundary.modes)
    coefficients = unknowns[mapping.grid_node_count:]
    inner = np.where(boundary.bi == grid.nr-1, pr, pz)
    outer = np.where(boundary.bi == grid.nr-1, dr@coefficients, dz@coefficients)
    return float(np.linalg.norm(inner-outer)/np.linalg.norm(outer))


def static_measurements():
    packet = CurrentPacket()
    exact_inside = packet.reference(INTERIOR_PROBES)[0]
    exact_outside = packet.reference(EXTERIOR_PROBES)
    result = {"source": asdict(packet), "interior_probes_r_z_m": INTERIOR_PROBES.tolist(),
              "exterior_probes_r_z_m": EXTERIOR_PROBES.tolist(), "exact_interior_psi": exact_inside.tolist(),
              "exact_exterior_fields": exact_outside.tolist(), "quadrature_order": 32,
              "quadrature_32_vs_48_relative": float(np.linalg.norm(packet.reference(EXTERIOR_PROBES, 48)-exact_outside)/np.linalg.norm(exact_outside)),
              "spatial": [], "modal": []}
    for cells in (8, 16, 32, 64):
        mapping, unknowns, rhs = static_packet(cells)
        result["spatial"].append({"radial_fluid_cells": cells, "matching_extent_m": 2., "modes": 12,
                                  "fluid_mesh": asdict(mapping.fluid), "field_grid": asdict(mapping.field),
                                  "relative_interior_flux_L2": float(np.linalg.norm(nodal_probes(mapping, unknowns)-exact_inside)/np.linalg.norm(exact_inside)),
                                  "relative_normal_derivative_jump": normal_derivative_jump(mapping, unknowns),
                                  "relative_solve_residual": float(np.linalg.norm(mapping.k@unknowns-rhs)/np.linalg.norm(rhs))})
    for modes in (1, 2, 4, 6, 8, 12, 16):
        mapping, unknowns, _ = static_packet(32, modes=modes)
        field = np.asarray(mapping.open_boundary.outside_field(unknowns[mapping.grid_node_count:], EXTERIOR_PROBES[:, 0], EXTERIOR_PROBES[:, 1]))
        result["modal"].append({"radial_fluid_cells": 32, "matching_extent_m": 2., "modes": modes,
                                "relative_exterior_B_L2": float(np.linalg.norm(field[1:]-exact_outside[1:])/np.linalg.norm(exact_outside[1:])),
                                "relative_exterior_psi_L2": float(np.linalg.norm(field[0]-exact_outside[0])/np.linalg.norm(exact_outside[0]))})
    return result


def pulse_measurement(cells, extent=2, *, boundary="multipole", modes=12, cfl=.2):
    solver = coil_pulse(cells, extent=extent, cfl=cfl, vacuum_boundary=boundary, multipole_modes=modes)
    history = [solver.diagnostics()]
    for time in np.linspace(0, .1, 21)[1:]:
        solver.advance(float(time)); history.append(solver.diagnostics())
    row = {"radial_cells": cells, "extent": extent, "boundary": boundary, "modes": modes,
           "configuration": configuration(solver), "diagnostics": solver.diagnostics(), "history": history}
    return row


def all_measurements():
    result = {"static": static_measurements(), "pulse": []}
    print("Completed independent asymmetric field and modal refinement", flush=True)
    cases = [(8, 2, "multipole", 12, .2), (16, 2, "multipole", 12, .2),
             (24, 2, "multipole", 12, .2), (32, 2, "multipole", 12, .2),
             (16, 3, "multipole", 12, .2), (16, 4, "multipole", 12, .2),
             (16, 2, "multipole", 1, .2), (16, 2, "multipole", 4, .2),
             (16, 2, "multipole", 8, .2), (16, 2, "multipole", 20, .2),
             (16, 2, "multipole", 12, .1),
             (16, 2, "fixed", 12, .2), (16, 3, "fixed", 12, .2), (16, 4, "fixed", 12, .2)]
    for cells, extent, boundary, modes, cfl in cases:
        result["pulse"].append(pulse_measurement(cells, extent, boundary=boundary, modes=modes, cfl=cfl))
        print(f"Completed pulse: nr={cells}, extent={extent}, {boundary}, modes={modes}, CFL={cfl}", flush=True)
    return result
