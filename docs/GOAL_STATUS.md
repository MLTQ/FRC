# Active objective: a credible experimental FRC simulation

The goal is a numerically verified and physically justified simulation in which
we can vary geometry, plasma and prescribed coil waveforms to investigate FRC
formation. Passing the existing stationary-conductor tests is not completion.

## Completion evidence still required

1. Conservative axisymmetric compressible dynamics: mass, momentum, pressure/
   total energy, magnetic advection and resistivity, with compatible div B.
2. Verified plasma–vacuum coupling and suitable chamber/end/outer boundaries;
   no undisclosed conducting or massive numerical atmosphere.
3. A defensible transport/ionization model for the selected regime, with local
   scale/collisionality checks; Hall/two-fluid or kinetic treatment if required.
4. Independent tests of waves, shocks, magnetic forces, diffusion/heating and
   equilibrium preservation, including convergence and energy budgets.
5. Converged formation experiments with separate reversal, topology, force
   balance, lifetime and model-validity criteria. Robustness to time step, mesh,
   boundary and transport choices; no claim of stability from axisymmetry alone.
6. Scriptable and interactive experiments with complete reproducibility records
   and diagnostics that do not overstate the supported conclusions.

## Current evidence

The initial turn was progress: Stage 0/1 were implemented with 33 passing tests
and three waveforms plus a 27-case sweep. Their topology durations are not
converged and their illustrative argon parameters raise missing-physics warnings.

The first continuation implemented a conservative dynamic r-z MHD kernel that
passes independent cylindrical fast-wave, 2-D flux-advection, exact Sod shock,
resistive-decay/heating, axis-balance, divergence and temporal-order checks.
The next continuation added limited second-order reconstruction of fluid and
staggered magnetic traces, exact regular-axis moments, total-pressure source
integration and improved one-sided wall gradients. Smooth wave/advection errors
are much smaller; axis force accuracy improves, and a particular straight-field
polynomial equilibrium is held to roundoff. First-order comparison remains
selectable and is rerun on identical grids.
The complete suite now has **62 passing tests**. At identical resolution the
fast-wave error is 216 times smaller and the 2-D flux-advection error 129 times
smaller than the first-order reference on the reported benchmark grids.

An independent continuum Grad–Shafranov fixture now tests genuinely curved,
closed poloidal flux with reversed Bz. It starts in equilibrium, not from coil
formation. Its artificial motion decreases with refinement, but it is not held
exactly. A lower-pressure p0=0.1 Pa stress case previously failed reconstruction
on coarse meshes despite positive cell pressure. That failure is now fixed by
local admissibility reconstruction, with consistent geometric-source integration.
Cell averages and magnetic traces are unchanged by this limiter; there are no
cell floors or energy repairs. The original failure remains reproducible with
the explicit disable switch. Limiting affects 12.10% and 2.19% of accepted
cell-stages at nr=12/24, and zero at nr=48/96 in that fixture.
The default p0=1 Pa case completes across refinement after correcting the wall
derivative. Smooth temporal refinement verifies RK2; the quarter-period wave's
limiter-switching time ratios are separately retained and less regular.
For p0=1 Pa, the finest 192×384 equilibrium run has maximum artificial speed
0.00095046 times the reference Alfvén speed after 0.1 s, with flux relative L2
drift 2.1181e-5. The initial force volume RMS converges at approximately order
1.5, so this remains a measured limitation rather than exact equilibrium balance.
Fixed-grid time refinement barely changes that drift, confirming that spatial
error dominates this particular case. Mass/energy and face divergence remain
at floating-point accuracy. Generated figures were visually inspected.

Additional equilibrium stress cases down to p0=1e-5 Pa finish conservatively,
but their coarse-grid artificial motion is not physical accuracy evidence.
An independently solved magnetized double rarefaction remains positive and
conservative with unchanged uniform B. Its bulk L1 errors converge, while the
central pressure is still 114% too high at 640 axial cells. The admissibility
limiter is inactive there: this exposes underlying flux-discretization error,
not a pressure floor. It is tracked as an accuracy requirement, not hidden by
the passing positivity checks. No full-update positivity theorem is claimed.

Measured evidence is in `docs/MHD_VERIFICATION.md`,
`docs/EQUILIBRIUM_VERIFICATION.md`, `docs/POSITIVITY_VERIFICATION.md`, and their
JSON/plot outputs. The limiter derivation is `docs/POSITIVITY.md`. These remain
bounded verification solvers, not the coupled coil-driven formation experiment.

## Current development step

Next, connect dynamics to the independent electromagnetic/vacuum subsystem
with a measured source-work budget. The bounded reference core and its
verification issue are complete in their declared scope; remaining cold-flow
accuracy, curved-equilibrium stress balance and stronger full-update checks
are tracked separately and must be addressed before formation claims. Transport
applicability can be developed independently. Do not
simply enable `Experiment(stage=2)` with periodic/metal-wall boundaries.

One physically relevant coupling route is a pre-ionized gas filling the actual
chamber, with fluid only inside the chamber and a true magnetic vacuum outside.
That avoids adding an artificial massive exterior. It still requires validated
quartz/end boundary treatment and must address rarefaction, changing conductivity
and formation of an internal plasma–vacuum interface if a detached torus forms.
It is an implementation candidate, not a completed or validated model.

## Completion audit at this checkpoint

| Requirement | Authoritative evidence | Status |
|---|---|---|
| Conservative dynamic core | mhd modules, independent tests, MHD/positivity reports | Implemented and verified for stated bounded reference; positive traces without cell floors; full-update robustness and cold-flow accuracy remain bounded by tests |
| Physical source/plasma/vacuum coupling | Original Stage 2 guard still raises; beads frc-vacuum | Incomplete |
| Regime-appropriate transport / Hall physics | Fixed eta/sigma only; beads frc-transport | Incomplete |
| Full independent numerical verification | Waves, shocks, heating, advection, J×B, curved equilibrium, admissibility and double rarefaction data | Partial; measured curved-equilibrium and cold-pressure error; no formation validation |
| Converged FRC formation | Only older stationary-conductor results, unconverged topology duration | Not demonstrated |
| Reproducible experimental interface | Stage 0/1 CLI/playback; dynamic benchmark scripts and metrics | Partial; no coupled dynamic experiment UI |

Beads issues `frc-goal`, `frc-vacuum`, `frc-transport`, `frc-accuracy`, and
`frc-formation` preserve the remaining scope; `frc-mhd` records the completed
bounded reference core milestone. No external process is being
waited on. The optional request for target plasma parameters is not a blocker
to numerical verification or development.

The full objective remains active until all completion evidence is present.
