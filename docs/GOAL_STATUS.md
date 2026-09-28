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
The bounded-reference checkpoint had **62 passing tests**. At identical resolution the
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

## Dynamic vacuum coupling checkpoint

The complete suite now has **74 passing tests** (18.08 s at this checkpoint).
`VacuumCoupledMHD` advances compressible gas inside a fixed, gas-filled chamber
with reflecting slip walls/ends and magnetically transparent boundaries. A
larger magnetic-only exterior satisfies the elliptic induced-current constraint;
there are no exterior density, pressure or conductivity degrees of freedom.
An adjoint magnetic force/advection pair, implicit midpoint and gas Euler Strang
splitting account for induced self-field energy, gas energy, prescribed-source
work and positive physical Ohmic heat. This response-energy budget excludes
external-field self/mutual energy and is not a coil-circuit energy estimate.
Failure leaves accepted state, time and cumulative work intact. A source time
cap is explicit because fluid CFL cannot discover arbitrary short pulses.

Independent compact-field J cross B errors decrease from 0.2707 to 0.006437 in
volume L2 as radial resolution increases 8→64, including axis cells. Passive
advection errors decrease 0.008325→0.002155→0.0005438 at 12/24/48 cells; Gaussian
resistive-diffusion errors decrease 0.007174→0.001822→0.0004573. Fixed-grid time
halving reduces momentum error by 4.010 and 4.050. A separate infinite-space
circular-loop integral distinguishes mesh error from magnetic-box error:
at 32 cells/m, zero-flux outer extents 2/4/8 m give 5.429%, 0.6791%, 0.08840%
probe error. With independent exact boundary data the 32-cell/m error is 0.03344%.

Two prescribed coils now induce current, gas acceleration and heat in this
reference. Across the reported dynamic cases the maximum relative response-
energy residual is 1.63e-14; mass change is at most 1.42e-16. Compatible face
divergence is near roundoff. **Pulse flow is not yet converged:** final peak speed
increases 0.0661→0.1473→0.1798→0.1999 m/s at 8/16/24/32 radial cells. The magnetic
boundary also matters: at 16 cells, extent 2→4 m changes final source work by
about 19%, and 3→4 m still changes it by about 3.8%. A small solver residual is
therefore not physical boundary accuracy. Figures were inspected; full inputs,
histories, source fingerprint and versions are retained with the report.

See `docs/COUPLED_NUMERICS.md`, `docs/COUPLED_VERIFICATION.md` and
`outputs/coupled_verification/metrics.json`. Earlier reports remain historical
measurements of their recorded source, not reruns with this different solver.

## Current development step

Prioritize a verified open magnetic boundary or demonstrated domain convergence,
resolution of source-driven skin layers and stronger compression/magnetic-shock
checks for the split method. Energy pairing alone does not prove local momentum
or shock accuracy. Replace the fixed gas-filled interface with a defensible
evolving plasma–vacuum treatment before interpreting detached torus formation.
Validate actual chamber/end conditions, transport and Hall applicability. The
fixed-chamber reference is useful measured progress, not completion of
`frc-vacuum`. `Experiment(stage=2)` stays gated.
Cold-flow pressure accuracy, curved-equilibrium stress balance and stronger
full-update robustness remain tracked separately. Transport applicability can
be developed alongside numerical coupling without selecting arbitrary closure.

## Completion audit at this checkpoint

| Requirement | Authoritative evidence | Status |
|---|---|---|
| Conservative dynamic core | mhd modules, independent tests, MHD/positivity reports | Implemented and verified for stated bounded reference; positive traces without cell floors; full-update robustness and cold-flow accuracy remain bounded by tests |
| Physical source/plasma/vacuum coupling | Fixed-chamber split reference and COUPLED_VERIFICATION; original Stage 2 guard still raises | Partial: smooth/work checks pass; finite-domain sensitivity, free interface, end conditions and magnetic shocks remain |
| Regime-appropriate transport / Hall physics | Fixed eta/sigma only; beads frc-transport | Incomplete |
| Full independent numerical verification | Waves, gas shocks, heating, advection, J×B, curved equilibrium, admissibility, rarefaction and new vacuum/coil budget checks | Partial; measured cold-pressure/equilibrium errors and coupled pulse/boundary sensitivity; no formation validation |
| Converged FRC formation | Only older stationary-conductor results, unconverged topology duration | Not demonstrated |
| Reproducible experimental interface | Stage 0/1 CLI/playback; dynamic and coupled benchmark scripts, histories and metrics | Partial; no coupled dynamic experiment UI |

Beads issues `frc-goal`, `frc-vacuum`, `frc-transport`, `frc-accuracy`, and
`frc-formation` preserve the remaining scope; `frc-mhd` records the completed
bounded reference core milestone. No external process is being
waited on. The optional request for target plasma parameters is not a blocker
to numerical verification or development.

The full objective remains active until all completion evidence is present.
