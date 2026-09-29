"""Reproduce the archived conservative study with its explicit global limiter."""
from pathlib import Path
import hashlib
import json
from .conservative_measurements import measure, provenance


def protocol():
    # The tracked study records are the reproducible input manifest. Read all
    # inputs before overwriting any requested output; never infer them from a
    # current solver default which may have changed since the reference study.
    archive = json.loads(Path("outputs/conservative_verification/metrics.json").read_text())
    return [(r["name"], r["kind"], r["radial_cells"], {**r["input_options"], "flux_limiter": "global"})
            for r in archive["cases"]]


def all_measurements(directory, *, reuse=False, names=None):
    cases = protocol()
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    source = provenance()
    source["reference_driver_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    source["reference_protocol_sha256"] = hashlib.sha256(json.dumps(cases, sort_keys=True).encode()).hexdigest()
    if names is not None and set(names)-{c[0] for c in cases}:
        raise ValueError("Unknown historical reference case")
    rows = []
    for name, kind, cells, options in cases:
        if names is not None and name not in names:
            continue
        path = directory/f"{name}.json"
        cached = json.loads(path.read_text()) if reuse and path.exists() else None
        if cached and cached.get("provenance") == source:
            row = cached
            print(f"Reused global reference {name}", flush=True)
        else:
            print(f"Starting global reference {name}", flush=True)
            row = measure(kind, cells, **options)
            row.update(name=name, provenance=source, input_options=options)
            path.write_text(json.dumps(row, indent=2, allow_nan=False)+"\n")
        rows.append(row)
    return {"provenance": source, "cases": rows}
