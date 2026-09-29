"""Separate gas-end and magnetic-matching controls for the driven boundary study."""
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
import hashlib
import json
import numpy as np
from plasma_bench.geometry import Grid
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.electromagnetics import VacuumField
from plasma_bench.coupled import CoupledState, VacuumCoupledMHD
from plasma_bench.coupled_hydro import gas_from_primitive
from plasma_bench.coupled_benchmarks import coil_pulse
from plasma_bench.conservative_benchmarks import from_split_initial
from .boundary_order_probe import install_probe
from .boundary_pulse_measurements import snapshot, source_record
from .verification_records import write_record


def measure(method, half_length, cells=64, directory="outputs/boundary_verification/ends"):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    name = f"{method}_length{half_length}_{cells}"
    start = perf_counter()
    # The reference supplies exactly the same source waveforms and plasma inputs.
    reference = coil_pulse(4, vacuum_boundary="multipole")
    mesh = MHDMesh(1., float(half_length), cells, 2*half_length*cells)
    grid = Grid(3., 3., 3*cells+1, 6*cells+1)
    source = VacuumField(grid, reference.sources.coils, bias=reference.sources.bias)
    split = VacuumCoupledMHD(mesh, source, CoupledState(gas_from_primitive(mesh, 1., 1.), np.zeros(grid.shape)), reference.parameters)
    solver = from_split_initial(split, energy_quadrature="staggered", boundary_poynting="variational", time_order=3)
    if method == "strip":
        install_probe(solver)
    provenance = source_record()
    provenance["end_controls_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    row = {"name": name, "method": method, "radial_cells": cells, "completed": False,
           "provenance": provenance, "parameters": asdict(solver.parameters),
           "fluid_mesh": asdict(mesh), "field_grid": asdict(grid),
           "source_configuration": {"bias_t": source.bias, "coils": [asdict(c) for c in source.coils]},
           "boundary_order_probe": getattr(solver, "boundary_order_probe", None),
           "setup_seconds": perf_counter()-start, "history": [solver.diagnostics()]}
    start = perf_counter(); fields = {}; times = []
    partial = directory/f"{name}.progress.json"
    try:
        for i, end in enumerate(np.linspace(0., .1, 21)[1:], 1):
            solver.advance(float(end)); row["history"].append(solver.diagnostics())
            row.update(diagnostics=solver.diagnostics(), advance_seconds=perf_counter()-start)
            if i in (8, 12, 20):
                fields.update({f"t{i}_{k}": v for k, v in snapshot(solver).items()})
                times.append({"prefix": f"t{i}", "time_s": solver.time})
            write_record(partial, row)
            print(name, solver.time, solver.steps, flush=True)
        row["completed"] = True
    except (ValueError, RuntimeError) as error:
        row["failure"] = str(error)
    row.update(diagnostics=solver.diagnostics(), advance_seconds=perf_counter()-start)
    path = directory/f"{name}.npz"; np.savez_compressed(path, **fields)
    row["fields"] = {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                     "times": times, "units": {"primitive": "rho kg/m3, vr/vz m/s, p Pa, Br/Bz T",
                     "total_flux": "Wb/(2pi)", "current": "A/m2", "numerical_eta_over_physical": "dimensionless"}}
    write_record(directory/f"{name}.json", row); partial.unlink(missing_ok=True)
    return row


def run(directory="outputs/boundary_verification/ends", *, reuse=False):
    rows = []; source = source_record()
    source["end_controls_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for method in ("standard", "strip"):
        for length in (1, 2):
            path = Path(directory)/f"{method}_length{length}_64.json"
            row = json.loads(path.read_text()) if reuse and path.exists() else None
            if not row or not row["completed"] or row["provenance"] != source:
                row = measure(method, length, directory=directory)
            rows.append(row)
    return rows


if __name__ == "__main__":
    run()
