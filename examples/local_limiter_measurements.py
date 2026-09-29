"""Persist local/global limiter controls without overwriting historical studies."""
from pathlib import Path
import argparse
import hashlib
import json
from .conservative_measurements import measure, provenance


def cases():
    result = [(f"shock_local_{n}", "shock", n, {}) for n in (24, 48, 96, 192, 384)]
    result += [(f"shock_global_{n}", "shock", n, {"flux_limiter": "global"}) for n in (96, 192)]
    result += [("shock_half_cfl", "shock", 192, {"cfl": .1}),
               ("shock_long_chamber", "shock", 96, {"half_length": 2}),
               ("shock_no_limiter", "shock", 96, {"positivity_flux_limiting": False})]
    result += [(f"pulse_{n}", "pulse", n, {}) for n in (16, 32, 64)]
    result += [("pulse_half_cfl", "pulse", 32, {"cfl": .1}),
               ("pulse_far_matching", "pulse", 16, {"extent": 4})]
    result += [(f"{kind}_{n}", kind, n, {}) for kind in ("advection", "diffusion") for n in (12, 24, 48)]
    legacy = {"energy_quadrature": "cell", "boundary_poynting": "riemann", "time_order": 2}
    return [(name, kind, cells, {**legacy, **options}) for name, kind, cells, options in result]


def source_record():
    record = provenance()
    record["local_measurements_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return record


def run(directory="outputs/local_limiter_verification/cases", *, reuse=False, names=None):
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source = source_record()
    if names is not None and set(names)-{c[0] for c in cases()}:
        raise ValueError("Unknown measurement case")
    rows = []
    for name, kind, cells, options in cases():
        if names is not None and name not in names:
            continue
        path = directory/f"{name}.json"
        cached = json.loads(path.read_text()) if reuse and path.exists() else None
        if cached and cached.get("provenance") == source:
            row = cached
            print(f"Reused {name}", flush=True)
        else:
            print(f"Starting {name}", flush=True)
            row = measure(kind, cells, **options)
            row.update(name=name, provenance=source, input_options=options)
            path.write_text(json.dumps(row, indent=2, allow_nan=False)+"\n")
            d = row["diagnostics"]
            print(f"{name}: completed={row['completed']}, t={d['time_s']:.6g}, "
                  f"response error={d['response_energy_balance_relative']:.6g}, "
                  f"max affected fraction={d['maximum_limiter_cell_fraction']:.6g}", flush=True)
        rows.append(row)
    return {"provenance": source, "cases": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse", action="store_true")
    parser.add_argument("--case", action="append")
    args = parser.parse_args()
    run(reuse=args.reuse, names=args.case)
