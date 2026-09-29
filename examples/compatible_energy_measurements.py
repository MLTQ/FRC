"""Persist compatible magnetic-energy, vacuum-work and time controls."""
from pathlib import Path
import argparse
import hashlib
import json
from .conservative_measurements import measure, provenance
from .verification_records import write_record


def cases():
    result = [(f"shock_{n}", "shock", n, {}) for n in (24, 48, 96, 192, 384)]
    result += [("shock_rk2_96", "shock", 96, {"time_order": 2}),
               ("shock_rk2_192", "shock", 192, {"time_order": 2}),
               ("shock_rk2_half_cfl", "shock", 192, {"time_order": 2, "cfl": .1}),
               ("shock_half_cfl", "shock", 192, {"cfl": .1}),
               ("shock_long_chamber", "shock", 96, {"half_length": 2}),
               ("shock_cell_energy", "shock", 96, {"energy_quadrature": "cell"}),
               ("shock_riemann_boundary", "shock", 96, {"boundary_poynting": "riemann"}),
               ("shock_no_limiter", "shock", 96, {"positivity_flux_limiting": False})]
    result += [(f"pulse_{n}", "pulse", n, {}) for n in (16, 32, 64)]
    result += [("pulse_half_cfl", "pulse", 32, {"cfl": .1}),
               ("pulse_far_matching", "pulse", 16, {"extent": 4}),
               ("pulse_rk2", "pulse", 32, {"time_order": 2})]
    result += [(f"{kind}_{n}", kind, n, {}) for kind in ("advection", "diffusion") for n in (12, 24, 48)]
    compatible = {"energy_quadrature": "staggered", "boundary_poynting": "variational", "time_order": 3}
    return [(name, kind, n, {**compatible, **options}) for name, kind, n, options in result]


def source_record():
    source = provenance()
    source["compatible_measurements_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return source


def run(directory="outputs/compatible_energy_verification/cases", *, reuse=False, names=None):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source = source_record()
    if names is not None and set(names)-{c[0] for c in cases()}:
        raise ValueError("Unknown compatible-energy case")
    rows = []
    for name, kind, cells, options in cases():
        if names is not None and name not in names:
            continue
        path = directory/f"{name}.json"
        progress_path = directory/f"{name}.progress.json"
        cached = json.loads(path.read_text()) if reuse and path.exists() else None
        if cached and cached.get("provenance") == source and cached.get("input_options") == options:
            row = cached
            print(f"Reused {name}", flush=True)
        else:
            print(f"Starting {name}", flush=True)
            def progress(partial):
                partial.update(name=name, provenance=source, input_options=options)
                write_record(progress_path, partial)
                d = partial["diagnostics"]
                print(f"{name}: accepted t={d['time_s']:.8g}, steps={d['steps']}, "
                      f"response residual={d['response_energy_balance_residual_j']:.8g} J", flush=True)
            row = measure(kind, cells, progress=progress, **options)
            row.update(name=name, provenance=source, input_options=options)
            write_record(path, row)
            progress_path.unlink(missing_ok=True)
            print(f"Finished {name}: completed={row['completed']}, "
                  f"failure={row.get('failure')}", flush=True)
        rows.append(row)
    return {"provenance": source, "cases": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse", action="store_true")
    parser.add_argument("--case", action="append")
    args = parser.parse_args()
    run(reuse=args.reuse, names=args.case)
