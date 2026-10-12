# Full-duration signed trajectory refinement

## First refined trajectory outcome (2026-10-12)

The first F8/64,B32/64 negative run is terminal: it accepted 17 steps through
28.333us, then rejected step 18 with a Newton GMRES convergence failure.
Its resource guards passed; session36840 exited0 through the driver's recorded
failure path and has been collected. This is not a successful100us trajectory.
The independent density-only prefix observation and its complete fresh replay
are described in TRAJECTORY_REFINEMENT_VERIFICATION.md.

This outcome changes execution order: first reproduce the complete rejected
protocol using the new one-job rejection_replay_plan.json. The prepared
seven-job plan remains preserved and unexecuted. Do not launch it concurrently
or automatically duplicate its first replay. The serial controller itself is
now used/frozen, with live status in outputs/trajectory_rejection_campaign.
After the replay is terminal and collected, inspect the reproduced failure
before selecting the opposite-sign, same-grid arithmetic and time controls.
All full-duration and full physical/conservation gates below remain required.

Continue from c912fdc, verified progress: all251 checkpoint evidence files,
168 frozen sources and the user's poster are unchanged. All preceding workers
are terminal/collected. The full FRC objective and all12 acceptance beads stay
open. No unit tests; complete physical workflows and independent audits only.

The radial16/N128 signed-rate amplitude campaign is complete in its stated
scope. Resume the separate pending radial1/N64 contact evolution benchmark.
Use the unchanged split_mode_evolution full nonlinear conservative fluid,
separate species entropy, finite-mass current and magnetic equations. The
independent continuum eigenmode initializes the signed perturbation and
supplies exp(lambda*t) audit targets only; no later state, forcing, growth fit
or projection is substituted into evolution. The benchmark growth is
8763.74925436694/s, distinct from the high-radial-mode rate study.

Keep the filled periodic argon fixture, fixed impermeable radial wall, constant
bias and physical masses. The original wall constraints, Gauss stages,
endpoint/induction residuals, positivity and nonlinear budgets remain. Absolute
stage tolerance1e-12 and relative1e-10 are unchanged. No explicit numerical
diffusion, resistivity, filtering, floors, atmosphere or energy repair is added.
The existing spatial approximation and Gauss method retain their unresolved
dispersion/diffusion and near-grid growth errors; acceptance is not accuracy.

Run both epsilon +/-1e-4 trajectories through100us on F8/64,B32/64, with60 and
120 steps. Every actual accepted stage/endpoint must be saved. A failed or
interrupted trajectory retains its last accepted state and cannot pass the
full-duration gate. Each new complete protocol requires a fresh full replay,
including all raw state/stage bytes and non-runtime metrics. Use fresh paths;
freeze new control code before its first execution. No prior source is edited.

First run the finer60-step negative trajectory and its fresh replay, then the
positive trajectory/replay and full independent opposite-sign audit. Actual
outcomes determine subsequent diagnosis; do not silently skip failed cases or
shorten duration. The120-step pair and replays remain required for time control.

The preserved F4/64,B16/64,60-step pair uses older induction arithmetic:
outputs/matched_mode_guarded_ne4_z64_t60_minus and _plus. It cannot be called
a bitwise implementation comparison to the split inverse. Run same-grid split
controls and full replays before attributing any cross-grid change to mesh
refinement. Independently audit both saved pairs using the same actual split
RHS for electric fields, retaining original source/arithmetic differences.
The earlier4.31611% final-current error and near-grid growing response remain
open; the refined grid is not presumed accurate.

Independent q16/q20 audits must reconstruct all fields, both pressures, current,
full electric field, growth projection, original continuous energy, particles,
species entropy and positivity. Normalize energy drift by the independent
initial perturbation-energy scale as well as background energy. Compare full
physical differences and mesh/time sensitivity, not growth projection alone.
Audit/comparison workflows require complete fresh replays. Local charge,
enthalpy and Poynting budgets are still separate unfinished gates.

## Resource budget

Keep8GiB own cumulative RSS and256MiB free disk guards, including final extrema,
and only one large process at a time. Source-derived retained modal arrays
are0.0579GiB and60-step secant history0.0397GiB for F8/64,B32/64; these exclude
other matrices, nonlinear/audit workspaces and Python, so are not peak bounds.

At start, free disk is410.156MiB. A preserved F4/64,B16/64,60-step trajectory
occupies14.5254MiB of compressed raw files. Doubling fluid/field radial counts
suggests roughly30MiB per new60-step trajectory, a planning estimate only.
The first primary plus full replay fits estimated headroom; measure actual
outputs before queuing later cases and account for Git storage. The entire
60/120-step campaign is not presumed to fit current disk headroom. Verified
byte-identical completed replays may share immutable storage, retaining every
path/hash. Never delete prior science or user files to make a run appear viable.

The first accepted F8/64,B32/64 step stores536064 uncompressed array bytes and
498195 compressed bytes; initialization stores180736 uncompressed bytes.
The60-step payload estimate is31.1909MiB. The prepared serial controller reserves
1.1 times uncompressed payload plus8MiB,42.3100MiB here, above the256MiB guard.
This is a planning margin, not a proven maximum; it excludes Git objects.

Originally prepared, not executed as a seven-job plan:
outputs/trajectory_refinement_diagnostics/t60_campaign_plan.json. After the
current manual primary is terminal and its actual process exit is collected,
the plan queues its full replay, the refined positive primary/replay and both
F4/64,B16/64 split-arithmetic primaries/replays. It stops before an unbudgeted
next job or after a resource/process error. Terminal numerical failures retain
their separately declared replays and never count as full-duration successes.
The120-step pair remains separate work after measured storage is reassessed.
The observed first failure supersedes that original launch order as stated above.

This is an evolution verification fixture, not formation or a lab prescription.
Refined-grid amplitude noise and1e-11 rejection remain open. Heated moving exact
vacuum, physical transport/ionization/model validity, real prescribed coils and
ends, converged formation and the interactive experiment UI remain required.
Neither this benchmark nor axisymmetric closed flux establishes3-D stability
or particle confinement. Keep research and issue export local.
