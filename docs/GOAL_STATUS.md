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

That coupling checkpoint had **74 passing tests** (18.08 s).
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

The magnetic boundary requirement has progressed with the verified infinite
exterior described below. The new shock study exposes false magnetic reversal
in the centered split and measures an explicit vanishing diffusion option.
An unsplit conservative fluid/UCT implementation now couples through the same
vacuum and preserves the tested interior frozen-flux relation. Prioritize
localizing its full-update flux limiter and extending shock/interface-work
verification, alongside source-driven skin-layer resolution. Its global limiter
is active on finer cold-boundary grids and affects otherwise smooth interior
accuracy. The historical split diffusion option remains a comparison. Replace the
fixed gas-filled interface with a defensible
evolving plasma–vacuum treatment before interpreting detached torus formation.
Validate actual chamber/end conditions, transport and Hall applicability. The
fixed-chamber reference is useful measured progress, not completion of
`frc-vacuum`. `Experiment(stage=2)` stays gated.
Cold-flow pressure accuracy, curved-equilibrium stress balance and stronger
full-update robustness remain tracked separately. Transport applicability can
be developed alongside numerical coupling without selecting arbitrary closure.

## Infinite magnetic exterior checkpoint

That checkpoint had **83 passing tests** (20.90 s). The preceding turn
made verified progress, with no external job remaining live at this checkpoint.
`vacuum_boundary="multipole"` now represents the induced magnetic field to
infinity using regular-axis decaying Legendre modes. Boundary nodal traces,
interior magnetic energy and an independently integrated exterior energy form
one symmetric matrix. Exterior nodes and modal amplitudes satisfy stationary
vacuum rows; modes have no fluid mass, current density or resistivity. The
midpoint response-energy budget includes energy beyond the field grid.
The historical fixed zero-flux boundary remains an explicit comparison.

Independent displaced toroidal-current quadrature excites both odd and even
modes. Interior flux relative L2 errors at 8/16/32/64 fluid radial cells are
0.01257, 0.0008744, 0.0001845 and 0.00004552. Interior/exterior normal-derivative
mismatch decreases 0.01653→0.003882→0.0009413→0.0002318, including matching-surface
corners. Exterior B error at fixed 32 cells drops from 32.3% with one mode to
about 0.0364% with twelve; higher modes then reach the spatial-error plateau.
A separate analytic dipole and independent all-space volume integral check
normalization, axis fields and the sign/magnitude of exterior magnetic energy.

All fixed/open pulse comparisons were rerun from current code. At 16 fluid radial
cells and 12 multipoles, moving the matching rectangle from 2 to 4 m changes
source work by **0.00206%**. The corresponding fixed-wall change is **23.646%**,
using the 4 m result as denominator for both percentages. The prior roughly 19%
figure used the 2 m fixed-wall result as denominator. Open pulse energy residuals
are at most 2.21e-14 relative to initial response energy; mass change is at most
2.83e-16 and modal stationarity residual at most 6.32e-13. Beyond-grid magnetic
energy is positive and included. No exterior fluid or state repair is introduced.

The pulse still is not spatially converged: open-boundary final peak speeds are
0.05919, 0.13227, 0.16262 and 0.18119 m/s at 8/16/24/32 cells. At 16 cells,
halving CFL changes source work by about 0.046%, larger than the measured boundary
surface effect. More modes alone cannot fix fluid or time-discretization error.
The current single-center expansion requires a matching rectangle enclosing the
entire chamber's source sphere; efficiency for long chambers remains to evaluate.

Derivation, measured results, inspected figure, complete inputs/histories and
source fingerprint are in `docs/OPEN_BOUNDARY.md`,
`docs/OPEN_BOUNDARY_VERIFICATION.md`, and
`outputs/open_boundary_verification/metrics.json`. This verifies the open
**magnetic** exterior, not a moving free plasma surface, gas outflow, transport,
magnetic shocks or FRC formation. `frc-vacuum` and the full objective stay active.

## Magnetic shock and implicit solver checkpoint

That checkpoint had **99 passing tests** (25.17 s). That continuation made
verified progress; all measurement/report jobs completed, with none left live.

An independent strong perpendicular shock now has exact Rankine–Hugoniot jumps
and steady cylindrical branches. Finite-difference momentum and invariant checks
validate the reference independently of the numerical kernels. Both actual 2-D
solvers are tested in an interior annulus before remote boundaries dominate.
The centered coupled method creates a **false Bz reversal** at nr=48, reaching
Bz/B1=-0.027288 even with a relative response-energy residual around 1e-14.
Such magnetic structure must not be interpreted as FRC topology.

The explicit `shock_diffusion` option adds a nonnegative compressive diffusivity
that scales as h² in smooth flow and h at shocks. It is off by default and
separately records numerical magnetic heat, without changing physical Ohmic
heat or putting fluid/resistivity in vacuum. No cell floors or energy repair are
introduced. With C=1, magnetic-profile L1/B1 errors at nr=24/48/96/192 are
0.330402, 0.218542, 0.106434, 0.0581267; density errors are 0.299483, 0.187281,
0.0960794, 0.046949 kg/m³; pressure errors are 1.93626, 1.06167, 0.574753,
0.190947 Pa. No false sign reversal occurs in those C=1 windows. However,
frozen B/rho error is **not monotonic**, and nr=192 still overshoots Bz/B1 to
3.6293 versus the physical shock value 3.53234. The conservative bounded
reference preserves that special B/rho invariant to roundoff. More diffusion
is not uniformly more accurate; the coefficient study retains that tradeoff.

At nr=48, moving the inner nonanalytic extension changes pressure/B profiles by
0.3302%/0.08888% in relative L2. Doubling chamber length changes profiles below
7.1e-9; a shorter 0.01 s wall-taper control gives density change 4.86e-5.
CFL halving changes the B profile by 1.43e-5, then 3.45e-6. These are measured
controls at that grid/window/time, not general free-boundary verification.
Maximum coupled response-energy residual is 1.94e-14, mass change 2.19e-16,
and face divergence 1.25e-16 T/m. Numerical magnetic heat remains explicitly
reported (0.433 J at nr=192 for the whole chamber); physical eta and heat are
zero. Entropy-producing shock dissipation need not vanish with grid spacing.

The midpoint solve now defaults to a cached-vacuum Picard iteration at tolerance
1e-13. It solves the same implicit equations and retains atomic failure. The
original direct solve remains selectable. Direct/Picard tests cover both
magnetic boundaries and independent physical/numerical diffusion choices. Their
nr=24 shock profiles agree within 9.9e-13 relative L2; measured advance time is
6.96 times faster for Picard in that fixture. This is a CPU wall-time measurement,
not a GPU or general real-time claim.

The inspected six-panel figure, 20 completed cases, complete configuration,
initial/final diagnostics, profiles, environment and source fingerprints are in
`docs/SHOCK_VERIFICATION.md`, `docs/SHOCK_NUMERICS.md`, and
`outputs/shock_verification/metrics.json`. Old reports retain historical hashes.
Next, investigate the split method's inconsistent magnetic/fluid numerical flux
and local momentum/weak-shock accuracy, while preserving the verified vacuum
and work identities. Continue the unresolved cold rarefaction, equilibrium,
skin-layer, free-interface, transport and formation requirements below.

## Conservative fluid/CT coupling checkpoint

The current full suite has **109 passing tests** (32.92 s). This continuation is
verified progress: all 22 measurement jobs completed, including 20 successful
runs and two deliberately disabled-method failure controls. No job remains live.

`ConservativeVacuumMHD` advances mass, momentum, total fluid-cell MHD energy and
upwind constrained induction together, with shared Rusanov/UCT diffusion. It
retains the exact discrete induced vacuum and infinite magnetic exterior, with
no exterior fluid or conductivity. Vacuum-side magnetic traction and Poynting
flux include interface-current work. Positive numerical diffusivity D=a h/4 is
assigned only to plasma boundary nodes and uses the symmetric stiffness current;
its diagnosed heat is separate from physical eta J². This is declared numerical
interface stabilization, not transport. Physical and numerical electric fields
enter all adjacent total-energy fluxes, including transport along the boundary
layer, with no double-counted Joule source. An independent manufactured Poynting
identity checks that a prescribed electric curl in initially uniform B changes
magnetic energy without a spurious instantaneous thermal-energy derivative.

The original strong shock fixture reaches 0.02 s at nr=24/48/96/192. Interior
frozen B/rho defects are 2.85e-6, 9.22e-8, 4.45e-13 and 9.27e-14, with no false
field reversal. B-profile L1/B1 errors are 0.309278, 0.203389, 0.104325 and
0.0746278. Cold boundary cells require blending complete conservative/CT fluxes
toward an admissible first-order Euler proposal: 0/88, 0/176, 6/372 and 154/810
Euler stages are limited. Minimum high-order weight is 0.000332 at nr=192.
No cell pressure/density floor or energy repair is used. A failed low-order
proposal still raises atomically; no unrestricted positivity theorem is claimed.
The disabled limiter fails at 0.019111 s on nr=96; disabling interface diffusion
and limiting fails at 0.0100 s on nr=24. These failures remain in the data.

The invariant improvement is not universal profile improvement: at nr=192 the
new B error 0.07463 exceeds the historical split C=1 error 0.05813, and pressure
error is 0.36880 Pa versus 0.19095 Pa. Global limiting caused by boundary cells
reduces accuracy in the whole domain, so localization is the immediate numerical
task. Shock time halving at nr=96 changes B/density profiles by 0.02976% and
pressure by 0.07051%; doubling chamber length at nr=48 changes profiles below
2.70e-8. The original initial data and independent reference were retained.

Fluid-cell energy flux accounting is at roundoff, but independently calculated
gas-plus-full-space-induced-field response energy is not forced exact. Relative
shock response errors decrease 1.996e-4, 8.248e-5, 4.050e-5, 2.873e-5. This
explicit finite-grid error is distinct from the historical split work identity.
Mass changes are at most 2.19e-16 and compatible divergence 2.09e-16 T/m across
completed cases. All quantities remain independently diagnosed.

The same prescribed-coil pulse runs at nr=8/16/32/64. Final response error
normalized to **actual source work** decreases in magnitude 2.722%, 1.108%,
0.3406%, 0.1021%. At nr=64, physical Ohmic heat is 0.082397 J and numerical
interface heat is 0.0014725 J (about 1.8% of physical heat). Source work still
changes 0.125487→0.121817 J from nr=32→64, and peak speed increases
0.06981→0.12737→0.17880→0.20744 m/s over the four grids. The driven flow remains
unconverged. No pulse stages require flux limiting. CFL halving at nr=32 changes
work by 0.001741%, and magnetic matching extent 2→4 at nr=16 changes work by
0.002050%. Removing numerical interface diffusion at nr=32 changes work by
0.6537% and peak speed by 1.4146%, so its effects must not be ignored.

Independent smooth advection errors at nr=12/24/48 are 0.05585, 0.009702 and
0.001778; Gaussian diffusion errors are 0.007897, 0.002115 and 0.0004833, with
no flux limiting. Harmonic-polynomial tests independently verify vacuum-side
interface traces and their radial convergence. Tests cover no vacuum numerical
resistance, rest/bias preservation, cold-boundary atomic failure, local energy
conservation and independent source-budget refinement. The inspected figure,
full inputs/histories/profiles/failure records and fingerprints are in
`docs/CONSERVATIVE_COUPLING.md`, `docs/CONSERVATIVE_VERIFICATION.md` and
`outputs/conservative_verification/metrics.json`.

This remains a fixed gas-filled chamber reference. Free plasma motion and
actual end conditions, transport/ionization/Hall validity, converged formation
and a dynamic experiment UI remain unfinished. Original Stage 2 stays gated.

## Completion audit at this checkpoint

| Requirement | Authoritative evidence | Status |
|---|---|---|
| Conservative dynamic core | mhd modules, independent tests, MHD/positivity reports | Implemented and verified for stated bounded reference; positive traces without cell floors; full-update robustness and cold-flow accuracy remain bounded by tests |
| Physical source/plasma/vacuum coupling | Fixed-chamber split and conservative/CT references; COUPLED, OPEN_BOUNDARY and CONSERVATIVE verification reports; original Stage 2 guard still raises | Partial: magnetic/work/interface refinement measured; global flux limiter, source-layer resolution, moving free interface and actual end conditions remain |
| Regime-appropriate transport / Hall physics | Fixed eta/sigma only; beads frc-transport | Incomplete |
| Full independent numerical verification | Waves, gas/perpendicular magnetic shocks, heating, advection, J×B, equilibrium, trace/full-update admissibility, rarefaction, interface traces, source budgets and infinite-exterior checks | Partial; conservative/CT removes tested frozen-flux defect, but global limiting and interface diffusion affect accuracy; cold-pressure/equilibrium, broader shocks and pulse-resolution requirements remain |
| Converged FRC formation | Only older stationary-conductor results, unconverged topology duration | Not demonstrated |
| Reproducible experimental interface | Stage 0/1 CLI/playback; dynamic and coupled benchmark scripts, histories and metrics | Partial; no coupled dynamic experiment UI |

Beads issues `frc-goal`, `frc-vacuum`, `frc-transport`, `frc-accuracy`, and
`frc-formation` preserve the remaining scope; `frc-mhd` records the completed
bounded reference core milestone. No external process is being
waited on. The optional request for target plasma parameters is not a blocker
to numerical verification or development.

The full objective remains active until all completion evidence is present.
