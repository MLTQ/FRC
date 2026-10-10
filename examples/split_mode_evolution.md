# Complete end-to-end split-inverse trajectory

Runs the frozen matched_mode_evolution accepted-step protocol with guarded
Gauss and SplitInverseWallPlasma. Only induction arithmetic changes. The full
actual nonlinear state, physical masses, stage/endpoint/induction bounds,
positivity, duration and saving/failure semantics remain unchanged. Source
identity explicitly includes this opt-in dispatch and its dependencies.

Uses the established first-radial-mode N64 contact benchmark, distinct from
the current high-wavenumber radial16,N128 spatial study. A successful full
trajectory here does not close that separate spatial gate. The independent
reference supplies only initialization and audit targets, never later states.

All temporary serial module bindings are restored on exit, including failure.
Use a separate process for parallel variants. Verify through full signed
trajectories, physical/conservation audits and complete stage/state replay;
the user's end-to-end-only preference excludes new or existing unit tests.
No prior frozen driver changes or silently adopts this arithmetic.
