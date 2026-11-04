# hops_campaign.py

First HOPS (**Harmonically Oscillating Plasma Sphere**) actual 3D kinetic
experiment, plus source-matched DC and fixed-period controls. Does not select
an optimum or demonstrate fusion. Reuses the positive-cage solver; no FRC
solver, running job or immutable campaign source is edited.

`main` freezes source with the existing campaign helper and executes its worker
from that copy. Each case has additional module provenance and full states.
The campaign also copies its input mesh before any child starts.
Only one electrostatic replay runs at a time, with one CPU/BLAS/Numba thread.

First complete +3 kV source-off 20 ns startups at 10 and 5 ps, same initial
sample/mesh and a 2 ns voltage response. Record particle-throughput-normalized
energy residual, charge closure and losses. Refuse to start long runs if the
5 ps maximum residual exceeds 5% or substantially worsens. This explicit
development gate does not establish spatial-force/phase accuracy.

Then actual 0–5 µs at 5 ps, 600 macros/species, n0=1e11 m^-3, D+ reference mass,
18/1 mm mesh, 1 mA / 2 eV ideal wall electron source, +3 kV ceiling, 20 ns
voltage response, 200 ns common on warmup. HOPS switches off at inward 0.8 R0,
on at outward R0, with a 10 ns radius sensor and 20 ns minimum dwell. DC stays
on. Fixed pulses use a declared 1 µs period and 50% duty, not a found resonance.
Switch intervals, full 3D anisotropy, live counts, coherent/residual energy,
well potential and raw source/absorption/field/particle ledgers stay visible.

CLI: `.venv/bin/python -m examples.hops_campaign outputs/hops_v1 --mesh
outputs/iec_mesh_dev/cage.npz`. Duration is in actual microseconds; timestep
remains picoseconds for real electron mass. No mass scaling, fixed electron
background, imposed well, symmetry filter, force smoothing or energy repair.
The historical air discharge is not modeled. Collisions, ionization, local
wire/phase accuracy and particle/mesh convergence remain open.

Campaign `STOP` prevents later cases. An active case's own `STOP` stops it at a
step boundary and retains its partial state. A failed case stops the queue.
Long cases are under `outputs/lab_runs/` so the existing app lists live replays.
