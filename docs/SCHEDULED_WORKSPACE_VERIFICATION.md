# Component stress enables separate field refinements

2026-10-12 continuation from `a67d240`. All workers are terminal and collected.
No unit tests were added or run. The full original FRC sandbox goal and all
12 acceptance beads remain open. These are complete coupled instantaneous
signed-rate calculations in a filled periodic cylinder, with an impermeable
radial wall and prescribed bias. They are not trajectories or formation runs.

## Verified workspace changes

`ScheduledCurrentTransport` moves independent velocity-gradient construction
after stress evaluation and release. Every other statement is AST-identical;
this changes call scheduling, unlike the preceding deletion-only source check.
Complete coarse F4/64, B16/64 and fine F64/256, B256/256 workflows reproduce every
original detached-cache input, rate, face array, strict solver diagnostic,
physical fingerprint and q16/20 norm bitwise. Both have fresh full simulation
replays and independently replayed original/new diagnostic comparisons.

Scheduling alone peaks at 6.481/6.476 GiB in primary/replay, compared with the
preceding detached variant's 6.446/7.056 GiB. This does not establish a robust
memory improvement. The required fine replay also records 25 read-only stage
observations. First transport raises cumulative peak from 4.894 to 5.766 GiB.
Second transport reaches 6.476 GiB before stress_rate entry, during full stress
construction. Volume and face audits add no new peak. These nested high-water
intervals are not exclusive or live allocation measurements. Every scientific
result still replays exactly; read the physical and observation archives together.

`ComponentDriftStress` evaluates `((m_reduced/e**2)*J_i)*J_j/n` only when the
unchanged stress adjoint requests a component. It preserves multiplication and
division association and evaluates off-diagonal entries separately, without
assuming floating-point symmetry. Only the stress-construction assignment
changes in the copied pair method; every other statement AST matches.
No actual density, physical mass, gradient, cylindrical connection, face
traction, work adjoint, quadrature, boundary or residual check changes.
No diffusion, filtering, floor, atmosphere, artificial mass or repair is added.
Default geometry and time integrators remain unchanged.

The component variant's complete coarse/fine workflows, original/new
comparisons and all fresh replays are bitwise identical. Its fine primary/replay
peak at 5.342/5.516 GiB. This measured headroom enables the guarded increments
below; it is not a guaranteed bound for another grid or process.

| Complete workflow | Primary peak (GiB) | Full replay peak (GiB) |
|---|---:|---:|
| Scheduling only, B256/256 | 6.481 | 6.476 |
| Component stress, B256/256 | 5.342 | 5.516 |
| Axial refinement, B256/320 | 6.901 | 7.162 |
| Radial refinement, B320/256 | 7.161 | 7.037 |

## Independent directional accuracy

All cases keep F64/256, the independent radial16/N128 reference, epsilon 1e-4,
physical masses and every original equation and acceptance bound. Only the
indicated field direction changes. Each new case completes both signed rates,
q16/20 volume and face audits, a fresh full simulation replay, and independently
replayed current-component and actual magnetic-rate curl reconstruction.
Every saved array and non-runtime metric replays bitwise.

Errors below are percentages of each component's own continuum L2 norm at q20:

| Physical rate | B256/256 | B256/288 | B256/320 | B320/256 |
|---|---:|---:|---:|---:|
| Radial current | 5.96939 | 3.05317 | 1.76199 | 5.90438 |
| Azimuthal current | 0.0190466 | 0.0164420 | 0.0157584 | 0.0121202 |
| Axial current | 1.18700 | 1.11402 | 1.09567 | 0.734849 |
| Radial magnetic field | 0.0979624 | 0.0491827 | 0.0266691 | 0.0976233 |
| Toroidal magnetic field | 0.547509 | 0.363014 | 0.302238 | 0.474179 |
| Axial magnetic field | 0.00994610 | 0.00518330 | 0.00311496 | 0.00976620 |

Axial refinement mainly improves radial current and magnetic response; radial
refinement improves axial/azimuthal current and radial velocity. The latter's
error drops from 0.183143% to 0.0983390%. Total-current error drops from
0.0878993% to 0.0552027%, still hiding weaker component errors.

Independent strong-curl radial-current error is 1.73644% for B256/320, versus
5.94247% at baseline. Weak/strong radial difference increases from 0.110758%
to 0.131541%, so the actual magnetic response improves without repairing the
weak Ampere identity. For B320/256, strong-curl axial-current error is 0.599842%
versus 1.00320% at baseline; its weak/strong difference drops from 0.345296% to
0.210047%. Strong curl is inside polynomial cells and excludes face distributions.
Initial face-representation diagnostics are retained separately from rate errors;
they do not establish evolved local charge/enthalpy/Poynting budgets.

The new axial/radial maximum residuals are 6.82535e-11 and 8.96429e-11 under
1e-10. Their largest instantaneous power defect divided by absolute exchanged
powers is 4.343e-9; this is not a local or time-integrated conservation claim.
New component q16/20 norm differences are at most 2.64e-16 absolute; physical
rate/axis/wall norms differ by at most 1.39e-15. New curl decompositions close
within 2.30e-16 relatively; reference squared norms agree across grids within
6.67e-16. Axial radial-current log slopes are 5.692 and 5.218 on the modest
256/288/320 increments. These are observations, not an established asymptotic
order or a resolution of small-amplitude cancellation sensitivity.

## Preserved interruption and resource bounds

The first axial replay was disk-guard interrupted during the plus-rate:
262144000 free bytes crossed the 268435456-byte guard. Own peak 7251804160
bytes stayed below 8 GiB. Worker 12516 exited 130 and was collected. There was
no recorded numerical rejection and no complete physical replay for that attempt.
Its partial inputs/rates, source and control files remain unchanged, with hashes
in `outputs/scheduled_workspace_diagnostics/axial320_replay_interruption.json`.
The fresh standalone replay succeeds from scratch. All successful final resource
observations satisfy both unchanged guards; host-wide disk changes have no
established single-process cause. No unrelated process or data was removed.

A verified storage inventory checks 1776 completed replay file paths against
unchanged verification summaries and raw SHA256. Atomic sharing of 852 additional
byte-identical inodes records 456679287 nominal duplicate bytes, preserving every
path and scientific byte; all paths are rehashed afterward. Observed free space
rises from 512000000 to 966656000 bytes. Earlier smaller sharing manifests also
remain. Every shared path is immutable: never overwrite it in place.

## Reproducibility and next action

Twelve replay groups comprise six complete signed simulation replays, four
original/new comparison replays and two independent curl-audit replays; eight
raw-file comparisons pass. Eight new modules have companions. All 49 new source
archives, 397 preceding archives and 22 older active source closures verify.
Historical initial-attempt archives retain their own historical source; they
are not asserted to match current code. Final source union: 164 files, SHA
`a7022b5b4d1a8f5b598c431e71f6865d51c178ca61eed43a33ef230cfaeeabf5`.

Metrics, manifests and checkpoint: `outputs/scheduled_workspace_diagnostics`.
Standalone PNG/PDF: `outputs/directional_workspace_report`; PNG visually checked.
All previous evidence and the user's poster remain unchanged. No publication or
remote is used. Follow NEXT_SESSION for refined amplitude/residual controls,
combined/further spatial refinement, full trajectories and local budgets.
Heated moving exact vacuum, transport/ionization and validity, prescribed real
coils/ends, converged formation and the experiment UI remain open. Provisional
lab inputs remain exploratory; no 3-D stability or confinement claim is made.
