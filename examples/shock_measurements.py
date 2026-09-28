"""Reproducible magnetic-shock profiles and independent solver/control comparisons."""
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import hashlib
import importlib.metadata
import json
import platform
import numpy as np
from plasma_bench.shock_benchmarks import standing_shock, initial_fields
from plasma_bench.coupled_hydro import gas_primitive
from plasma_bench.mhd_state import primitive


def provenance():
    digest = hashlib.sha256()
    for path in sorted(Path("plasma_bench").glob("*.py"))+[Path(__file__)]:
        digest.update(path.name.encode()); digest.update(path.read_bytes())
    return {"source_sha256": digest.hexdigest(), "python": platform.python_version(), "platform": platform.platform(),
            "versions": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "matplotlib")},
            "interpretation": "Local standing-shock verification; initialized shock, no FRC formation or physical transport claim",
            "window": {"r_m": [.44, .66], "z_m": [-.1, .1]},
            "error_measure": "Cylindrical overlap-volume weighted mean absolute cell-observable error; pressure recovered from independently averaged gas conserved variables",
            "time_units": "seconds", "physical_resistivity_ohm_m": 0.}


def measure(cells, *, end_time=.02, **options):
    start = perf_counter()
    solver, ref = standing_shock(cells, **options)
    setup_seconds = perf_counter()-start
    mesh = solver.mesh
    gas, flux = initial_fields(mesh, ref, inner=options.get("inner", .38), wall_taper=options.get("wall_taper", True))
    exact = gas_primitive(mesh, gas, ref.gamma)
    exact_bz = mesh.cell_fields(flux)[1]
    coupled = options.get("coupled", True)
    row = {"radial_cells": cells, "options": options, "requested_end_time_s": end_time,
           "reference": asdict(ref), "downstream_at_shock": list(ref.downstream),
           "mesh": asdict(mesh), "parameters": asdict(solver.parameters), "setup_seconds": setup_seconds,
           "boundary": "reflecting gas and infinite magnetic exterior" if coupled else "periodic axial and conducting radial wall",
           "field_grid": asdict(solver.mapping.field) if coupled else None,
           "initial_diagnostics": solver.diagnostics(), "completed": True}
    start = perf_counter()
    try:
        solver.advance(end_time)
    except (ValueError, RuntimeError) as error:
        row.update(completed=False, failure=str(error))
    row["advance_seconds"] = perf_counter()-start
    row["diagnostics"] = solver.diagnostics()
    if coupled:
        q = gas_primitive(mesh, solver.state.gas, ref.gamma)
        br, bz = mesh.cell_fields(solver.mapping.fluid_flux(solver.total_flux()))
    else:
        q = primitive(mesh, solver.state, ref.gamma)
        br, bz = q[..., 4], q[..., 5]
    radial = .5*np.maximum(np.minimum(mesh.rf[1:], .66)**2-np.maximum(mesh.rf[:-1], .44)**2, 0.)
    axial = np.maximum(np.minimum(mesh.zf[1:], .1)-np.maximum(mesh.zf[:-1], -.1), 0.)
    weight = radial[:, None]*axial[None, :]
    weight /= weight.sum()
    mean = lambda value: float(np.sum(value*weight))
    row["errors"] = {"density_L1_kg_m3": mean(abs(q[..., 0]-exact[..., 0])),
                     "radial_velocity_L1_m_s": mean(abs(q[..., 1]-exact[..., 1])),
                     "pressure_L1_pa": mean(abs(q[..., 3]-exact[..., 3])),
                     "Bz_L1_over_upstream_B": mean(abs(bz-exact_bz))/abs(ref.upstream_field),
                     "frozen_ratio_defect_L1": mean(abs(bz/ref.field_ratio-q[..., 0]))/ref.upstream_density,
                     "minimum_Bz_over_upstream_B": float((bz/ref.upstream_field)[weight > 0].min()),
                     "maximum_Bz_over_upstream_B": float((bz/ref.upstream_field)[weight > 0].max()),
                     "transverse_B_L1_over_upstream_B": mean(abs(br))/abs(ref.upstream_field)}
    average_z = lambda value: (value@axial/axial.sum()).tolist()
    row["profile"] = {"r_m": mesh.r.tolist(), "density": average_z(q[..., 0]), "radial_velocity": average_z(q[..., 1]),
                      "pressure": average_z(q[..., 3]), "Bz_over_upstream_B": average_z(bz/ref.upstream_field),
                      "exact_density": average_z(exact[..., 0]), "exact_pressure": average_z(exact[..., 3]),
                      "exact_Bz_over_upstream_B": average_z(exact_bz/ref.upstream_field)}
    return row


def all_measurements(directory, *, reuse=False):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source = provenance()
    cases = [(f"bounded_{n}", n, {"coupled": False}) for n in (24, 48, 96, 192)]
    cases += [(f"centered_{n}", n, {"shock_diffusion": 0.}) for n in (24, 48, 96)]
    cases += [(f"diffusive_{n}", n, {"shock_diffusion": 1.}) for n in (24, 48, 96, 192)]
    cases += [(f"coefficient_{c}", 48, {"shock_diffusion": c}) for c in (.5, 2.)]
    cases += [(f"time_{c}", 48, {"shock_diffusion": 1., "cfl": c}) for c in (.1, .05)]
    cases += [("inner_control", 48, {"shock_diffusion": 1., "inner": .37}),
              ("ends_control", 48, {"shock_diffusion": 1., "half_length": 2}),
              ("wall_tapered_short", 48, {"shock_diffusion": 1., "end_time": .01}),
              ("wall_untapered_short", 48, {"shock_diffusion": 1., "wall_taper": False, "end_time": .01}),
              ("direct_comparison", 24, {"shock_diffusion": 1., "midpoint_method": "direct"})]
    rows = []
    for name, cells, options in cases:
        path = directory/f"{name}.json"
        cached = json.loads(path.read_text()) if reuse and path.exists() else None
        if cached and cached.get("provenance") == source:
            row = cached
            print(f"Reused verified source/configuration: {name}", flush=True)
        else:
            row = measure(cells, **options)
            row.update(name=name, provenance=source)
            path.write_text(json.dumps(row, indent=2, allow_nan=False)+"\n")
            print(f"Completed {name}: success={row['completed']}, {row['advance_seconds']:.2f}s; B error={row['errors']['Bz_L1_over_upstream_B']:.6g}", flush=True)
        rows.append(row)
    return {"provenance": source, "cases": rows}
