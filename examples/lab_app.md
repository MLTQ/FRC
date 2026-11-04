# lab_app.py

## Purpose
Local-only unified scientific app at `http://localhost:8792`, serving a buildless
UI, experiment definitions, saved FRC viewers and owned electrostatic launches.
No cloud connection, publication, arbitrary shell execution or FRC job control.

## Components and contracts
- `Lab`: one active owned IEC child, validated numerical config, new output/log
  per launch. Children freeze source at entry. A default matching mesh can be
  reused from local cache; nondefault resolution builds its own mesh.
- `Handler`: GET `/api/state`; token/origin-checked POST `/api/launch` and
  `/api/stop`. Only allowlisted UI/output file types are served, with path
  traversal and directory-list rejection. Bound to loopback only.
  Host is restricted to localhost/127.0.0.1 at the selected port, including GET,
  to prevent a remote-origin hostname from rebinding onto the local service.
- GET `/api/observables?run=ID&frame=N&magnetic=0|1` derives full-population
  moments and native electric-field lines read-only from that saved state.
  Optional B is explicitly a regularized plasma-only contribution, not evolved
  total B. `ObservableStore` bounds caches and serializes this work. No NPZ files
  or arbitrary paths are exposed. Existing loopback host policy remains.
- `main`: `.venv/bin/python -m examples.lab_app [--port 8792]`.

Owned-run cancellation writes STOP; it never sends signals to unrelated jobs.
Partial evidence is retained. Stopping the HTTP server leaves physical children
running. Existing FRC viewers are read-only adapters under `/data/formation_display`.
The adapter serves the FRC viewer's native `.f16` volume buffers unchanged.
Extending the registry does not imply all solvers share equations or validity.
