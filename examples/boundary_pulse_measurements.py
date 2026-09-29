"""Source-driven convergence of the experimental boundary strip, with saved fields."""
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import argparse
import hashlib
import json
import numpy as np
from scipy.constants import mu_0 as MU0
from plasma_bench.coupled_benchmarks import coil_pulse
from plasma_bench.conservative_benchmarks import from_split_initial
from plasma_bench.vacuum_trace import interface_resistivity
from .boundary_order_probe import install_probe
from .conservative_measurements import provenance
from .verification_records import write_record


def protocol():
    cases = [(f"{method}_{n}", method, n, .2, 2) for method in ("standard", "strip") for n in (16, 32, 64, 128)]
    cases += [(f"{method}_half_time", method, 64, .1, 2) for method in ("standard", "strip")]
    cases += [(f"{method}_far_matching", method, 32, .2, 4) for method in ("standard", "strip")]
    return cases


def source_record():
    source = provenance()
    for name in ("boundary_pulse_measurements", "boundary_order_probe"):
        source[name+"_sha256"] = hashlib.sha256(Path(f"examples/{name}.py").read_bytes()).hexdigest()
    return source


def snapshot(solver):
    """Keep actual conserved/CT state and independently interpretable current."""
    m = solver.mapping
    q = solver._primitive(solver.state, solver.time)
    ku = m.k@m.interior(solver.state.induced_psi)
    current = np.divide(ku, MU0*m.node_area, out=np.zeros(m.node_count), where=m.active)
    numerical_eta = interface_resistivity(m, q, solver.parameters.gamma)
    return {"conserved": solver.state.conserved.copy(),
            "primitive": q.copy(), "total_flux": m.fluid_flux(solver.total_flux()).copy(),
            "current": m.fluid_flux(m.full(current)),
            "numerical_eta_over_physical": m.fluid_flux(m.full(numerical_eta/solver.parameters.resistivity))}


def measure(name, method, cells, cfl, extent, directory):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source = source_record()
    start = perf_counter()
    split = coil_pulse(cells, extent=extent, cfl=cfl, vacuum_boundary="multipole")
    solver = from_split_initial(split, energy_quadrature="staggered", boundary_poynting="variational",
                                time_order=3, max_timestep=.002*cfl/.2)
    if method == "strip":
        install_probe(solver)
    elif method != "standard":
        raise ValueError("Unknown boundary method")
    row = {"name": name, "method": method, "radial_cells": cells, "completed": False,
           "provenance": source, "parameters": asdict(solver.parameters),
           "fluid_mesh": asdict(solver.mesh), "field_grid": asdict(solver.mapping.field),
           "source_configuration": {"bias_t": solver.sources.bias,
                                    "coils": [asdict(c) for c in solver.sources.coils]},
           "boundary_order_probe": getattr(solver, "boundary_order_probe", None),
           "setup_seconds": perf_counter()-start, "history": [solver.diagnostics()]}
    start = perf_counter(); last_progress = start
    partial = directory/f"{name}.progress.json"
    fields = {}; saved_times = []
    def progress(current, force=False):
        nonlocal last_progress
        now = perf_counter()
        if force or now-last_progress >= 30:
            row.update(diagnostics=current.diagnostics(), advance_seconds=now-start)
            write_record(partial, row)
            print(name, "t=", current.time, "steps=", current.steps, flush=True)
            last_progress = now
    progress(solver, True)
    try:
        for i, end in enumerate(np.linspace(0., .1, 21)[1:], 1):
            solver.advance(float(end), on_step=progress)
            row["history"].append(solver.diagnostics())
            if i in (8, 12, 20):
                fields.update({f"t{i}_{k}": v for k, v in snapshot(solver).items()})
                saved_times.append({"prefix": f"t{i}", "time_s": solver.time})
            progress(solver, True)
        row["completed"] = True
    except (ValueError, RuntimeError) as error:
        row["failure"] = str(error)
    row.update(diagnostics=solver.diagnostics(), advance_seconds=perf_counter()-start)
    field_path = directory/f"{name}.npz"
    np.savez_compressed(field_path, **fields)
    row["fields"] = {"file": field_path.name, "sha256": hashlib.sha256(field_path.read_bytes()).hexdigest(),
                     "times": saved_times, "units": {"primitive": "rho kg/m3, vr/vz m/s, p Pa, Br/Bz T",
                     "total_flux": "Wb/(2pi)", "current": "A/m2", "numerical_eta_over_physical": "dimensionless"}}
    write_record(directory/f"{name}.json", row)
    partial.unlink(missing_ok=True)
    print("Finished", name, "completed=", row["completed"], flush=True)
    return row


def run(directory="outputs/boundary_verification/pulses", *, names=None, reuse=False):
    directory = Path(directory)
    source = source_record(); rows = []
    known = {c[0] for c in protocol()}
    if names is not None and set(names)-known:
        raise ValueError("Unknown pulse case")
    for args in protocol():
        if names is not None and args[0] not in names:
            continue
        path = directory/f"{args[0]}.json"
        row = json.loads(path.read_text()) if reuse and path.exists() else None
        if row and row["provenance"] == source and row["completed"]:
            field_path = directory/row["fields"]["file"]
            if not field_path.exists() or hashlib.sha256(field_path.read_bytes()).hexdigest() != row["fields"]["sha256"]:
                raise ValueError("Saved pulse fields do not match the recorded hash")
        else:
            row = measure(*args, directory)
        rows.append(row)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append")
    parser.add_argument("--reuse", action="store_true")
    args = parser.parse_args()
    run(names=args.case, reuse=args.reuse)
