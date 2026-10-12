# Resource-guarded complete nonlinear trajectories

`run(destination, control, spec, epsilon, steps, original=None)` executes the
unchanged split_mode_evolution protocol for100us at the original absolute
stage tolerance1e-12 and relative tolerance1e-10. The caller supplies separate
fluid/field meshes, signed nonzero epsilon and positive step count. The fixed
radial1/N64 continuum fixture supplies initialization and audit targets only.
All actual nonlinear fluid, species, current and magnetic states evolve; no
later state or forcing is taken from the reference.

The physical source identity remains that of the frozen evolution driver.
A separate control archive explicitly seeds the dynamic wrapper and actual
matched evolution, fixture, solver and replay source closures. No old Python
file is edited. Every original stage, endpoint, induction and positivity bound,
accepted state/stage file, progress update and terminal failure remains.

An own-process monitor keeps8GiB cumulative peak RSS and256MiB free disk limits.
It sends SIGINT on violation and exits only its own worker with130 if native
work has not yielded after five seconds. Recorded final extrema must also
pass.
Resource summaries are atomically refreshed every30 monitor observations;
they are status snapshots, not proof that a process is still live.
Resource interruption, terminal trajectory failure and successful full
duration are distinct. A failed run retains its last accepted step; it cannot
pass the full-duration gate. Ordinary solver failure is not retried at a
different step, tolerance or amplitude. Run one large process at a time.

`--original` requires a complete same-source protocol replay through the
existing all-state/stage verifier. Raw file hashes must also be identical, so
shape/value equality alone cannot hide dtype or signed-zero differences.
The accepted_full_trajectory_replay flag is false for a reproduced terminal
failure. Sharing storage is a separate action permitted only after a complete
verified replay. Fresh physical/control destinations are mandatory.

Full opposite-sign independent physical/conservation audits, same-grid
solver-arithmetic comparisons, mesh/time controls and their fresh complete
replays remain required. These filled periodic contact fixtures contain no
moving vacuum or physical ends and are not formation experiments. There is no
new diffusion, filter, floor, artificial mass, atmosphere or energy repair.
No unit tests. Resource estimates do not substitute for observed limits.
