"""Recorded shock, vacuum-interface and source-work accuracy of conservative MHD."""
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import hashlib
import importlib.metadata
import json
import platform
import numpy as np
from plasma_bench.conservative_benchmarks import from_split_initial
from plasma_bench.shock_benchmarks import standing_shock, initial_fields
from plasma_bench.coupled_benchmarks import coil_pulse, passive_advection, passive_diffusion, compact_flux, gaussian_diffusion
from plasma_bench.coupled_hydro import gas_primitive


def provenance():
    digest = hashlib.sha256()
    for path in sorted(Path("plasma_bench").glob("*.py"))+[Path(__file__)]:
        digest.update(path.name.encode()); digest.update(path.read_bytes())
    return {"source_sha256": digest.hexdigest(), "python": platform.python_version(), "platform": platform.platform(),
            "versions": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "matplotlib")},
            "interpretation": "Unsplit conservative fixed-chamber reference; independent response-energy error, no formation claim"}


def shock_errors(solver, reference):
    mesh = solver.mesh
    gas, flux = initial_fields(mesh, reference)
    exact = gas_primitive(mesh, gas, reference.gamma)
    exact_b = mesh.cell_fields(flux)[1]
    q = solver._primitive(solver.state, solver.time)
    radial = .5*np.maximum(np.minimum(mesh.rf[1:], .66)**2-np.maximum(mesh.rf[:-1], .44)**2, 0.)
    axial = np.maximum(np.minimum(mesh.zf[1:], .1)-np.maximum(mesh.zf[:-1], -.1), 0.)
    weight = radial[:, None]*axial[None, :]; weight /= weight.sum()
    mean = lambda value: float(np.sum(value*weight))
    errors = {"density_L1_kg_m3": mean(abs(q[..., 0]-exact[..., 0])),
              "pressure_L1_pa": mean(abs(q[..., 3]-exact[..., 3])),
              "Bz_L1_over_upstream_B": mean(abs(q[..., 5]-exact_b))/reference.upstream_field,
              "frozen_ratio_defect_L1": mean(abs(q[..., 5]/reference.field_ratio-q[..., 0]))/reference.upstream_density,
              "minimum_Bz_over_upstream_B": float((q[..., 5]/reference.upstream_field)[weight > 0].min()),
              "maximum_Bz_over_upstream_B": float((q[..., 5]/reference.upstream_field)[weight > 0].max())}
    average = lambda value: (value@axial/axial.sum()).tolist()
    profile = {"r_m": mesh.r.tolist(), "density": average(q[..., 0]), "pressure": average(q[..., 3]),
               "Bz_over_upstream_B": average(q[..., 5]/reference.upstream_field),
               "exact_density": average(exact[..., 0]), "exact_pressure": average(exact[..., 3]),
               "exact_Bz_over_upstream_B": average(exact_b/reference.upstream_field)}
    return errors, profile


def measure(kind, cells, *, cfl=.2, extent=2, half_length=1, **overrides):
    start = perf_counter()
    reference = None
    if kind == "shock":
        split, reference = standing_shock(cells, cfl=cfl, half_length=half_length)
        end, count = .02, 5
    elif kind == "pulse":
        split = coil_pulse(cells, extent=extent, cfl=cfl, vacuum_boundary="multipole")
        end, count = .1, 21
    elif kind in ("advection", "diffusion"):
        split = passive_advection(cells, cfl=cfl) if kind == "advection" else passive_diffusion(cells, cfl=cfl)
        end, count = (.05 if kind == "advection" else .03), 2
    else:
        raise ValueError(kind)
    solver = from_split_initial(split, **overrides)
    row = {"kind": kind, "radial_cells": cells, "parameters": asdict(solver.parameters),
           "fluid_mesh": asdict(solver.mesh), "field_grid": asdict(solver.mapping.field),
           "reference": asdict(reference) if reference else None, "completed": True,
           "requested_end_time_s": end, "setup_seconds": perf_counter()-start,
           "source_configuration": {"bias_t": solver.sources.bias,
                                    "coils": [asdict(coil) for coil in solver.sources.coils]},
           "initial_diagnostics": solver.diagnostics()}
    history = [solver.diagnostics()]
    start = perf_counter()
    try:
        for time in np.linspace(0, end, count)[1:]:
            solver.advance(float(time)); history.append(solver.diagnostics())
    except (ValueError, RuntimeError) as error:
        row.update(completed=False, failure=str(error))
        history.append(solver.diagnostics())
    row.update(advance_seconds=perf_counter()-start, history=history, diagnostics=solver.diagnostics())
    select = solver.limiter_cell_counts > 0
    i, j = np.nonzero(select)
    row["limiter_map"] = {
        "encoding": "Sparse cell indices (radial, axial); omitted cells have zero counts",
        "shape": [solver.mesh.nr, solver.mesh.nz], "indices": np.column_stack((i, j)).tolist(),
        "affected_stage_counts": solver.limiter_cell_counts[select].tolist(),
        "troubled_stage_counts": solver.limiter_troubled_counts[select].tolist(),
        "r_m": solver.mesh.r.tolist(), "z_m": solver.mesh.z.tolist(),
        "maximum_depth_from_gas_boundary_m": float(np.minimum(solver.mesh.radius-solver.mesh.r[i],
                                                               solver.mesh.half_length-abs(solver.mesh.z[j])).max()) if i.size else 0.}
    if reference:
        window = (solver.mesh.r[:, None] > .44)&(solver.mesh.r[:, None] < .66)&(abs(solver.mesh.z[None, :]) < .1)
        row["limiter_map"]["affected_reference_window_cell_stages"] = int(solver.limiter_cell_counts[window].sum())
    if reference:
        row["errors"], row["profile"] = shock_errors(solver, reference)
    elif kind in ("advection", "diffusion"):
        grid = solver.mapping.field
        r, z = grid.r[:, None], grid.z[None, :]
        exact = compact_flux(r, z, amplitude=1e-6, width=.5, shift=.025) if kind == "advection" else gaussian_diffusion(r, z, .03)
        select = (r < .7)&(abs(z) < .7)
        row["flux_relative_L2"] = float(np.linalg.norm((solver.state.induced_psi-exact)[select])/np.linalg.norm(exact[select]))
    return row


def all_measurements(directory, *, reuse=False):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source = provenance()
    cases = [(f"shock_{n}", "shock", n, {}) for n in (24, 48, 96, 192)]
    cases += [("shock_no_limiter", "shock", 96, {"positivity_flux_limiting": False}),
              ("shock_no_interface_diffusion", "shock", 24, {"interface_diffusion": False, "positivity_flux_limiting": False}),
              ("shock_first_order", "shock", 96, {"spatial_order": 1}),
              ("shock_half_cfl", "shock", 96, {"cfl": .1}),
              ("shock_long_chamber", "shock", 48, {"half_length": 2})]
    cases += [(f"pulse_{n}", "pulse", n, {}) for n in (8, 16, 32, 64)]
    cases += [("pulse_half_cfl", "pulse", 32, {"cfl": .1}),
              ("pulse_far_matching", "pulse", 16, {"extent": 4}),
              ("pulse_no_interface_diffusion", "pulse", 32, {"interface_diffusion": False})]
    cases += [(f"{kind}_{n}", kind, n, {}) for kind in ("advection", "diffusion") for n in (12, 24, 48)]
    rows = []
    for name, kind, cells, options in cases:
        path = directory/f"{name}.json"
        cached = json.loads(path.read_text()) if reuse and path.exists() else None
        if cached and cached.get("provenance") == source:
            row = cached
            print(f"Reused {name}", flush=True)
        else:
            row = measure(kind, cells, **options)
            row.update(name=name, provenance=source, input_options=options)
            path.write_text(json.dumps(row, indent=2, allow_nan=False)+"\n")
            print(f"{name}: completed={row['completed']}, t={row['diagnostics']['time_s']:.5g}, response error={row['diagnostics']['response_energy_balance_relative']:.4g}, limited stages={row['diagnostics']['limited_euler_stages']}", flush=True)
        rows.append(row)
    return {"provenance": source, "cases": rows}
