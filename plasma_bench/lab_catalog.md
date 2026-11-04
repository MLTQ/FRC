# lab_catalog.py

## Purpose
Shared local experiment registry and read-only adapter to existing FRC displays.
New experiments add registry entries and a solver adapter, not a new app shell.

## Components and contracts
- `catalog`: user-visible model names, explicit scope and physical launch fields.
- `validate`: allowlisted finite numeric inputs, bounded timestep/resource counts;
  rejects unknown fields and nonintegral particle counts. Bounds are local resource
  limits, not evidence that all values are physically/numerically resolved.
- `runs(root)`: reads lab run manifests, direct `iec_3d_*` experiments and existing formation `runs.json` without
  modifying/adopting FRC processes. URLs remain local to the server.

IEC defaults are exploratory, not an optimum. FRC launches stay in their existing
source-frozen campaign tooling; the app provides their saved visualization.
