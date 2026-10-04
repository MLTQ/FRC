# FRC virtual plasma bench — new-thread handoff

**Continuation update, 2026-10-04:** The user has resumed research in the new
chat. The regular species component is now implemented, with 416 passing tests,
24 bitwise-replayed cases and independent physical/local-enthalpy audits.
Read the latest GOAL_STATUS.md and NEXT_SESSION.md for the authoritative next
step: complete coupled moving trajectories and physical interface work under
frc-7mm.4. REGULAR_SPECIES_VERIFICATION.md records the new source and limits.
The user's provisional lab envelope is in RESEARCH_QUESTIONS.md. The original
handoff below is preserved as historical context; its "not implemented" species
status and instruction to start that component have been superseded.

---

Prepared 2026-10-04 at the user's request to stop this thread and continue in a
new one. **Do not continue research in the old thread.** The overall project
is incomplete; no completion claim is made. Resume only from the user's new
thread instruction. This document is a portable entry point, not a new task
scope or permission to publish anything.

## Objective and original scope

Build a numerically verified and physically credible simulation in which the
user can vary chamber/plasma/coil parameters and prescribed waveforms to
investigate formation of an FRC plasma torus. The broader product is a scriptable,
interactive plasma sandbox with progressively more capable physics models.
Preserve the full original brief in [ORIGINAL_PROJECT_REQUEST.md](ORIGINAL_PROJECT_REQUEST.md).
Initial axisymmetric stages are prerequisites; eventual 3-D/kinetic extensions
remain part of the longer-term roadmap, not something the present results prove.

Hardware stays abstract prescribed I(t)/B(t) waveforms. No pulsed-power supply
construction/design is requested. Use SI units; report numerical diffusion,
physical assumptions and model-validity limits explicitly. Field reversal,
closed flux, equilibrium, FRC-like transient, stability and confinement are
separate findings. A pretty render or a conserved grid energy is not sufficient.

## Repository and exact boundary

- Workspace: `/Users/max/Code/FRC`, branch `master`, no Git remote.
- Latest **research** checkpoint: `7b3c8b2` — Verify current transport and finish
  remote fluid audit. A later local commit contains this handoff only.
- At handoff inspection the only pre-existing uncommitted change was the user's
  `outputs/frc_showcase/poster.png`. Preserve it; it was excluded from our commits.
- The interrupted continuation performed **read-only inspection only**. No new
  species solver was implemented and no new numerical job was started after
  `7b3c8b2`. All research jobs launched by the preceding completed turn are terminal.
- Read `AGENTS.md` and module companion `.md` files before changing code. This
  project uses beads; its full local issue state is `.beads/issues.jsonl`.
- Beads executable: `/Users/max/.local/bin/bd --no-daemon`. Use `show <id>`;
  `sync --flush-only` saves locally. Do not invent a remote or publish outputs.
  Parent/child issue listings are not a reason to close the parent early.
- Existing beads/modular-docs skills apply. No sub-agent delegation is requested.

## Verified evidence to keep

1. **Current transport:** 16 degree/resolution/quadrature/density cases, all
   original-platform arrays/metrics replayed bitwise, with independent continuous
   reconstruction at quadrature orders 16 and 20. **408 tests passed in 228.47 s**.
   Source `2f682eece854d54b3b0d33026709f8d34617cc4633d8b3875a0f3e2b8e0d9540`.
   Read `CURRENT_TRANSPORT_NUMERICS.md` and `CURRENT_TRANSPORT_VERIFICATION.md`.
   Evidence: `outputs/current_transport_verification/` and
   `outputs/current_transport_diagnostics/` (logs, source checks, PNG/PDF report).
   The prior missing density-energy exchange is paired without adding heat.
   Local electric/force errors converge, including axis/wall regions.
2. **Limits remain measured:** grid total-power backward error <=3.04005e-14
   does not imply exact continuous energy conservation. Degree-three physical
   total-power error refines from 2.35289e-5 to 3.72857e-11 W. Strong continuity,
   div(J) and current jumps are not all identically zero. The prototype has no
   completed species transport, complete coupled evolution or formation result.
3. **Old remote fluid study:** all 14 trajectories complete and audited. The last
   nr32 contact completed 4116 steps; its initial state/rate and final accepted
   step replay bitwise, and order-18 physical integration agrees to 3.04540e-15.
   Evidence: `outputs/face_radial_remote_completion/` plus earlier local/remote
   records. `frc-5b5` is closed. The obsolete instruction to wait for/restart that
   run must not be followed.
4. **Species charge probe:** four resolutions at both quadrature orders replay
   bitwise. Source `c6a23e1abb844963d2eababc68a6b8c61a3d808ec8321dbee609b9bec2351c59`.
   Evidence: `outputs/species_charge_probe/`; read
   `REGULAR_SPECIES_COUPLING_DESIGN.md`. This is a diagnostic plus a proposed
   discrete construction, **not an implemented species integrator**.
5. All 13 earlier calculation source closures were checked unchanged. Archived
   sources, immutable raw data, failed attempts and diagnostic comparisons must
   remain reproducible. Use new modules/output directories for changed physics;
   never rerun over terminal evidence with different source.

## Remaining beads — all ten unfinished issues

| Bead | Status at handoff | Work still required |
|---|---|---|
| **frc-7mm.4** | in progress | **Start here.** Implement and verify shared species pressure/current-entropy transport, separate local temperatures, and enthalpy/Poynting work; then original-parent and unforced coupled tests. |
| **frc-7mm.2** | in progress | Complete regular-fluid/current/magnetic coupling and full moving evolution; component work identities alone do not finish this issue. |
| **frc-7mm.1** | in progress | Resolve/quantify remaining physical convective-energy, momentum and fluid order/spectral limitations; retain grid versus continuous-energy distinctions. |
| **frc-7mm** | in progress | Parent integration: conservative finite-mass moving fluid/thermal/current system, true total-energy flux and atomic failure, with space/time convergence. |
| **frc-accuracy** | in progress | Cold-flow pressure accuracy, curved-equilibrium stress balance and stronger low-pressure multidimensional robustness without floors or repairs. |
| **frc-vacuum** | in progress | Heated moving plasma/vacuum interface, external sources, physical chamber/end boundaries and source-work accounting; exactly zero exterior plasma. |
| **frc-transport** | in progress | Defensible collisions/conductivity, ionization/neutrals/radiation as required, and local Hall/two-fluid/kinetic applicability for the selected regime. |
| **frc-7mm.5** | open | Validate complete moving trajectories, then coordinate adoption, exact-vacuum, physical-closure and formation gates. The old remote run is already finished. |
| **frc-formation** | open | Converged prescribed-waveform formation experiments; distinguish reversal, topology, balance, lifetime and validity; verify mesh/time/boundary/transport robustness and reproducible controls. |
| **frc-goal** | in progress | Overarching acceptance: a credible experimental FRC sandbox, scriptable and interactive controls/sweeps, truthful diagnostics and complete reproducibility. Keep open until all requirements are proved. |

Recently completed: **frc-7mm.3** (current-transport audit) and **frc-5b5**
(old near-axis force verification and final trajectory audit). Do not reopen
these just because the larger integration remains unfinished.

## First implementation task in the new thread

Read the species design and these companions/sources:
`plasma_bench/regular_field_map`, `regular_current`, `current_transport`,
`physical_derivatives`, `regular_coordinate_fluid`, `compatible_geometry`,
`examples/species_charge_probe`, and `current_transport_sampling`.

Use the actual regular state n,q=n*vr/r,l=n*vtheta/r,mz,Ce,Ci; true angular
momentum is 2*x*l with x=r²/2. Transverse axis traces are dependent, not new
velocity degrees of freedom. Preserve the physical electron/ion masses and
common actual density/continuity. The pressure identity is

    grad(p_s)/n = grad(gamma_s*h_s) - h_s*grad(s_s),
    E_pressure = [alpha grad(p_i) - beta grad(p_e)]/(e n).

The saved candidate constructs pressure and entropy terms with compatible
potential gradients, whose current work vanishes by the curl/gradient identity.
That helps algebraic conservation but does not establish local consistency,
interface transport or stability. The probe explicitly measures the potential
interpolation error; no global energy subtraction is an acceptable substitute.

Implement in separate modules and verify against independent original pressure
and entropy rates with nonuniform density, independently varying temperatures
and discontinuous element traces. Measure global/axis/wall/max errors, both
species entropy integrals, charge/particle continuity, angular momentum, actual
physical energy and higher-quadrature convergence. Check frozen-current entropy
transport stability as well as work identities. Then combine components and
run complete unforced parent wave/contact trajectories, including time refinement
and enthalpy/Poynting boundaries. Do not replace verified sources prematurely.

Unimplemented note from the interrupted reasoning: a skew-adjoint version of
an entropy advection operator may be useful for preserving its frozen-density
quadratic invariant. This has **no code or measured evidence**. If explored,
derive the matching pressure term again (the compensating thermodynamic
potential changes); do not assume the existing nonskew proposal is already a
verified stable transport scheme.

## Execution, external state and tools

Typical local environment:

```sh
cd /Users/max/Code/FRC
PYTHONDONTWRITEBYTECODE=1 VECLIB_MAXIMUM_THREADS=1 \
  MPLCONFIGDIR=/private/tmp/plasma-bench-mpl MPLBACKEND=Agg \
  .venv/bin/python -m pytest -q
```

The test log already records the 408-pass baseline. Rerun relevant tests after
actual changes; do not spend time repeating long completed studies without a
reason. Existing study CLIs protect frozen directories from overwrite.

The user authorized SSH and data transfer to **m@192.168.0.202**. The completed
old run is `/home/m/Code/frc-tool-comparison/face-radial-20260930`. The original
PIDs were absent on October 4, and terminal evidence was audited/retrieved.
No pending simulation is expected from that study. Revalidate process identity
and terminal/checkpoint files before any future operation; an observation timeout
is never a reason to launch a duplicate. The work there used CPU NumPy/SciPy;
GPU acceleration has not been demonstrated. See ENGINE_COMPARISON.md for the
existing tooling investigation; COMSOL is not required by the current project.

The original request mentions an RTX 4090, but the new thread should inspect
actual remote hardware before selecting a GPU runtime. Target experimental
parameters and suitable papers/transport data remain inputs to physical model
selection; do not invent validated closures from illustrative argon fixtures.

## Existing visual artifacts

- `outputs/current_transport_diagnostics/convergence.png` and `.pdf`: latest
  verified transport/conservation report.
- `outputs/frc_showcase_checkpoint/index.html`, `poster.png`, `orbit.mp4`:
  interactive 3-D current and magnetic fields, based on frozen coupling case 07.
  The movie is **camera motion around one saved fixture**, not plasma evolution.
  Its displayed 402-test label is historical, not the latest solver status.
- Preview, if needed: `.venv/bin/python -m examples.showcase_checkpoint --serve --port 8779`.
  Check whether the port is already in use before starting another preview.
- Preserve the older initialized-FRC showcase and the user's modified poster.

## Suggested first message for the new thread

Continue the FRC project in /Users/max/Code/FRC. Read AGENTS.md,
docs/THREAD_HANDOFF.md and docs/GOAL_STATUS.md, then inspect bead frc-7mm.4.
Preserve the full original scope and verified archives. Implement and verify
the remaining species pressure/entropy coupling before complete moving-plasma
validation. Keep the overarching FRC goal open until its full evidence gates pass.
