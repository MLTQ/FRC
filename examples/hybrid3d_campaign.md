# hybrid3d_campaign.py

Runs the six user-requested physical initial-fill hybrid shots, with bounded
concurrency, immutable source copies and honest progress records. It changes
resolution and timestep, not equations, cutoff, symmetry or physical inputs.
No unit tests or verification hierarchy are added.

- `run_campaign(OUT)` uses 10×10×10 mm, eight ions per gas cell, 0.5 ns ion
  steps, coupled banks and conservative electron energy, from actual 0–40 µs.
  Full fields are saved every 0.1 µs, with diagnostics at the same interval.
  A startup gate now waits for baseline1.3µs, past the old1.141µs failure,
  before releasing a second worker. A baseline failure defers queued shots.
  At most two workers, two Numba threads each by default; BLAS uses one.
- `freeze_source` copies modules and viewer assets, records SHA256 hashes and
  the current git diff. Children import from that copy, using absolute data
  paths. Subsequent shared-checkout edits cannot change queued shot code.
- `progress` reads atomic physical histories, source work, energy residual,
  explicit cutoff/empty-node heat removal and loaded-ion counts. Accounting
  error above max(5 J, 1% of source work) invalidates that owned child; no
  physical state is rescaled. Other configurations may still run.
- Completed shots export to the shared viewer with distinct names and labels;
  coarse results remain. Failed exports retain the complete simulation.
  `DONE` means every queue entry is terminal, not that every run succeeded.
- Less than 10 GB disk free defers queued shots without deleting data.
  Interrupting the runner terminates only its own active children. It does
  not pause or modify other jobs. No automatic retry or checkpoint restart.
  SIGTERM is handled through the same owned-child cleanup.

Usage: `.venv/bin/python -m examples.hybrid3d_campaign outputs/hybrid3d_shots_r10`
from the project directory. Runs can take many hours; `campaign.json` and
per-shot logs are the current evidence. Numerical losses and grid/particle
sensitivity remain open in `docs/HYBRID_3D.md`; this is not a validated
confinement or lifetime prediction.

Actual 10 mm startup thread scaling is recorded in `docs/HYBRID_3D.md`:
two threads give 1.209× one-thread throughput; four add only 4% over two
in the measured 30-step interval. This supports keeping two threads per
current worker. Parallel configurations affect campaign throughput, whereas
threads affect a single shot; extra configurations can also contend for
memory bandwidth. The runner has no live worker/thread reconfiguration.
Source-frozen shots do not load later shared-checkout edits. New queues can
explicitly select the Metal backend below.
`--gpu 1` selects Metal transport/CFL evaluation; `--gpu 2` additionally
selects atomic PIC moments. The default remains CPU. GPU execution requires
GPU-capable process access, and the copied source and actual setup identify
the mode. Earlier source-frozen queues are not changed by these new flags.
The old r10 CPU queue was stopped after its positivity failure; corrected queues
must use new output names rather than overwriting that scientific evidence.
`--cases uniform_h2` can run a single matched CPU control while the main
queue covers all six; it does not adopt or change another queue's workers.
Current corrected queues: `outputs/hybrid3d_shots_r10_cfl_metal/` uses
`--gpu 2 --workers 2 --threads 2 --gate-us 1.3`, and
`outputs/hybrid3d_shots_r10_cfl_cpu_control/` uses
`--gpu 0 --workers 1 --threads 2 --cases uniform_h2 --gate-us 1.3`.
The latter validates precision against the same corrected energy operator;
neither full-shot convergence nor late-shot GPU throughput is yet established.
`--arrays` enables the float64 loop-fused electron stencils for a newly
created queue and records that choice in its coordinator, rows and physical
setup. It works with CPU and Metal. Default remains off; existing frozen
queues retain their current code and progress. One-shot CLI equivalent:
`arrays=1 thermal=1`. See `hybrid3d_array_check.md` for matched physical
replays against a pre-change source copy and measured timing limits.
