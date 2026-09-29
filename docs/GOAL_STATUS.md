# Active objective: a credible experimental FRC simulation

The goal is a numerically verified and physically justified simulation in which
we can vary geometry, plasma and prescribed coil waveforms to investigate FRC
formation. Passing the existing stationary-conductor tests is not completion.

Latest checkpoint: **257 passing tests; 63 completed compatible Hall cases**.
The prior axis error is corrected, and a separate compressible Hall solver now
evolves nonlinear flow/fields with implicit Hall/resistive induction and explicit
compatible magnetic energy. Physical B representation remains first order;
separate electron/ion thermal transport, finite-mass and vacuum/edge coupling
remain incomplete. See COMPATIBLE_HALL_VERIFICATION and TWO_TEMPERATURE_PLAN.

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

The separate fixed-material `MaterialCurrent` component now includes physical
electron inertia, spatially prescribed positive carrier density/nonnegative
resistivity, true magnetic vacuum and explicit electron drift kinetic energy.
Conditional Coulomb transport and local Hall/inertia/drift diagnostics are
integrated into reference experiment metadata. Prescribed conductivity and
moving dynamics are unchanged. The full suite has **192 passing tests (79.10 s)**.
All 35 current-component cases, an independent time audit, 43 conditional
states and 41 prescribed-edge samples are complete. Both final figures were
visually inspected; complete source/configuration/field hashes are checked.
No numerical job remains live at this checkpoint.

CURRENT_VERIFICATION.md separates continuum error, time order, physical current
energy and numerical false reversal. The unresolved de/R=.05, nr=8 screening
case produces negative Bz down to -0.0033032 of applied bias, despite a small
energy residual; the independent exact field is positive. Finer sampled fields
remove this artifact but nr=64 still has 7.3078% magnetic-field L2 error. The
inertial decay's nr=64 current/heat errors are 0.082887%/0.10486%; initial weak
surface-current error limits convergence and is retained, not hidden.

The original Ar+ reference has di/R=4.0972 and conditional electron mean-free-
path/R=0.11088. It does not justify ignoring Hall physics. Smaller physical
edge-gradient lengths raise electron inertia, carrier-drift and kinetic
requirements. The heated sharp-edge incompatibility in SURFACE_COMPATIBILITY
remains unresolved; the new linear fixed-density component does not fix it.

The next bounded coding step is the axisymmetric Hall/toroidal reference and
independent wave-dispersion verification in EXTENDED_PLASMA_PLAN.md. Its parent
species equations, current-energy force/transport terms, separate electron/ion
thermal evolution and justified edge treatment must precede complete moving
coupling. Spatial coil force/work checks remain in MOVING_COIL_PLAN.md. Actual
gas ends, mesh/shock robustness, converged formation and the experimental UI
remain required. The original low-pressure scope is retained; a denser/easier
fixture cannot silently replace it. `Experiment(stage=2)` stays gated and the
full objective remains active.

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

## Local complete-update limiter checkpoint

This continuation is verified numerical progress, not completion of the FRC
objective. The suite has **117 passing tests** (46.65 s). The default conservative
solver now chooses first-order shared face/corner fluxes only around inadmissible
Euler cells and rechecks newly troubled neighbors. Metric sources and actual
boundary energy fluxes follow that same local choice. Every flagged cell has
the independently checked low-order update, giving conditional finite
termination without a hidden global fallback. No floors or thermal repairs
were introduced; failed baseline proposals still raise atomically. The previous
global convex limiter remains an explicit, freshly reproduced control.

On the unchanged original shock, nr=192 B-profile L1 error decreases from
0.0746278 globally to 0.0593295 locally, and pressure L1 from 0.368802 to
0.334502 Pa. Frozen B/rho defect remains 8.06e-14. Only the two outermost radial
cell layers are directly affected; maximum stage support is 0.8545% of cells,
with zero direct replacements in the central independent-reference window.
Local limiting occurs more often (320 stages) but leaves distant fluxes intact.
The manufactured neighbor-chain test requires four successive detection passes;
the original nr=192 case needs at most two. Mass, fluid energy-flux accounting,
compatible divergence and vacuum stationarity remain near roundoff.

This improvement exposes a separate accuracy limitation. The nr=192 local
response-energy residual is +0.0140806 J (1.40031e-4 of initial response energy)
versus -0.00288898 J globally. The local residual is also **5.22% of the net
outward fluid energy**, a more demanding transfer scale, though that boundary
transfer is not external source work. Halving CFL raises the local residual to
0.0237089 J while central B/density profiles change only 0.001610% and pressure
0.002069%. Thus the interface evolution/energy is not time-converged simply
because the central shock is. Doubling chamber half-length at nr=96 changes
central profiles below 9.01e-9; the independent response budget remains separate.

The additional nr=384 run completes to the original 0.02 s end time: central
B/density L1 decreases to 0.0306517, pressure L1 to 0.175216 Pa, frozen-ratio
defect is 9.72e-13, and no reference-window cell is directly limited. Maximum
affected support is 0.4286% of cells, with at most two replacement passes.
However, response-energy error grows to **0.0575258 J, 30.20% of net outward
fluid energy**, and numerical interface heat rises from 0.06331 to 0.09848 J
between nr=192 and 384. The interface energy budget is therefore not converged.
Minimum density reaches 0.000388047 kg/m³ without a floor; 1575 accepted steps
take 2064 s after 101 s setup on the recorded CPU environment. The large-grid
run finished; no numerical process remains pending at this checkpoint.
All 21 case records (20 completed, one deliberately disabled-limiter failure)
and both energy audits have been source-fingerprint checked, and the final
scientific figure has been visually inspected.

The same coil pulses, time/matching controls, and smooth advection/diffusion
cases reproduce historical measured source work, response energy and analytic
errors exactly, with no local limiting. The earlier 0.1021% pulse error relative
to source work and unresolved peak-flow/skin-layer limitations remain.

An instantaneous response-rate identity now separates outward fluid power,
cell magnetic-energy rate and full-space nodal magnetic work. An independent
symmetric perturbation completes the magnetic vacuum and differences response
energy. It checks this diagnostic without changing any evolved state or budget.
The best symmetric-difference errors across the five snapshots stay below
2.62e-8 W at nr=96 and 3.65e-7 W at nr=192; probing leaves final diagnostics
exactly unchanged. At the final nr=192 state, the high-order defect rate is
-204.075 W, the first-order rate -2.83969 W, and the selected local rate
-7.23180 W. Thus this wall-state discrepancy is present in the spatial update,
well above diagnostic roundoff; its cumulative and timestep-sensitive behavior
still needs a consistent interface treatment.
The derivation is in `docs/INTERFACE_ENERGY_AUDIT.md`; measured rates, perturbation
errors, full controls, limiter maps and source fingerprints are retained in
`docs/LOCAL_LIMITER_VERIFICATION.md` and `outputs/local_limiter_verification/`.
The next interface task is to use that decomposition to improve consistent
vacuum/CT energy exchange without sacrificing local conservation or adding heat
repairs. Stronger multidimensional shocks, cold rarefaction, free surfaces,
actual ends, transport/ionization/Hall validity, formation and the dynamic
experiment UI all remain in scope.

## Compatible-energy and free-boundary-reference checkpoint

Positive staggered face-square magnetic energy and a natural variational
vacuum trace are now explicit solver options, with optional SSP RK3. Pressure,
initial energy, complete-update admissibility and diagnostics use the same
selected quadrature. Boundary work uses actual selected CT corner fields;
there are no thermal repairs or exterior fluid variables. Historical defaults
remain selectable and unchanged. Independent polynomial/force/work tests cover
unequal grid spacings, nonautonomous temporal order, third-stage atomic rejection
and mixed-boundary fallback. The pulse16 and advection12 legacy checks reproduce
historical source work, energy and analytic errors exactly.

The original shock completes nr=24/48/96/192/384. At nr=192 the response residual
falls from historical 0.0140806 J to 0.000596124 J; at nr=384 it falls from
0.0575258 J to **0.00173278 J**, or 0.728851% of net outward fluid transfer.
The fine B L1/B1 remains 0.0306581 and pressure L1 0.175312 Pa; frozen-ratio
defect is 7.07e-13. No central-reference-window cell is directly limited.
The finest run finishes 1408 steps with minimum density 0.000916625 kg/m³,
3423/4224 limited Euler stages, at most 0.48489% affected cells per stage and
three replacement passes. It takes 98.9 s setup and 2792.9 s advancement.
Its energy error still grows from nr=192, so full boundary convergence is
not established by this substantial improvement.

A per-step quadratic RK audit independently predicts the entire nr=192
response error to roundoff, with maximum stage spatial defect 2.20e-11 W and
maximum per-step identity error 1.13e-12 J. The constant-bias semidiscrete work
mismatch is resolved for this method; its finite-step error is not repaired.
At nr=24, CFL halving reduces RK3 energy error by factors 8.88 and 7.85.
The fine boundary is not time-converged: half-CFL gives 0.00100662 J with RK3
and 0.00218710 J with RK2. The full first-order spatial baseline instead gives
3.97637e-7 and 4.79246e-8 J, but degrades B L1 to 0.15656.

An explicitly experimental fixed first-order boundary strip retains the
second-order interior. Recorded coarse/fine controls test that hypothesis;
its forced order reduction is counted separately from any adaptive fallback.
It is not a production default. Source-driven skin layers, boundary accuracy
and finer shock meshes must be checked before adopting it.
The nr=64 unmodified compatible pulse residual is **0.00315436%** of applied
source work, versus historical 0.102119%. Peak flow remains mesh-sensitive.
All original pulse/smooth cases complete with inactive adaptive limiting.

The spherical FRC reference independently verifies true vacuum outside,
continuous magnetic traction with no sheet current, exact force balance and
full-space energy integrals. The existing topology detector finds one O-ring
and two axis nulls; separatrix-volume error falls from 0.219889% to 0.00374111%
over 24→192 intervals. This is initialized ideal equilibrium, not formation,
a density/transport closure, stability or confinement evidence. Its figure
has been visually inspected. See FREE_BOUNDARY_REFERENCE.md.

The full suite now has **135 passing tests (49.04 s)**. All 25 primary cases
(24 complete, one retained disabled-limiter failure), two first-order controls,
six boundary-strip probes and four RK audits are finished. Their source hashes
match the frozen engine/measurement fingerprint
`5247cc3d60224324a8f1f731b164f9fd41a16a8ea70a140af70ae1a302c022c2`.
Both scientific figures were visually inspected. No numerical jobs remain.
The nr=192 strip probe preserves B L1 0.059316 and gives response residual
3.97369e-7 J, decreasing to 4.78879e-8 J at half CFL; its 766 forced cells
and 1528 affected cells are explicitly reported separately from adaptive
fallback. Source-driven and finer-boundary checks remain open for this probe.
See COMPATIBLE_ENERGY.md and COMPATIBLE_ENERGY_VERIFICATION.md; full
source-matched records are in outputs/compatible_energy_verification.
The next numerical step is to verify the fixed boundary strip with driven
pulses and finer meshes, or improve the higher-order boundary closure without
losing conservation, CT and admissibility. Free-boundary evolution, actual
ends, transport/ionization/Hall validity, converged formation and the dynamic
experiment UI remain required. The Stage 2 formation gate stays closed.

## Boundary source, skin and end validation checkpoint

The current suite has **142 passing tests (57.97 s)**. A new independent
infinite-cylinder Bessel solution checks sinusoidal skin diffusion, regular-axis
and vacuum matching, Maxwell/Ohm identities and time-averaged surface/volume
Ohmic power. Fifteen recorded grid/time cases compare physical-only diffusion,
the production induction stencil and the fixed first-order boundary strip.
Both production choices converge in magnetic field; the surface-current
stencil has separate, slower pointwise accuracy. At nr=128 the maximum B L2
errors normalized by applied amplitude are 0.000934414 (standard) and
0.00150768 (strip), with about 7.22 cells per analytic skin depth.
This is stationary-conductor induction verification, not evolved gas dynamics.

All twelve full coil-pulse runs and four gas-end controls are complete, with
hashed lossless conserved/primitive/CT/current snapshots at 0.04, 0.06 and
0.10 s. The standard nr=16/32/64 budgets and peak speeds reproduce the previous
compatible study exactly. At nr=128, standard/strip final peak speeds are
0.221845/0.217821 m/s, but their maximum same-grid velocity L2 difference is
**9.76713%**. The 64→128 velocity changes are **7.05911% / 10.4227%**.
Thus the strip's favorable shock energy result has a material driven-flow
accuracy cost and it remains experimental; no production default is changed.
Density/pressure perturbations and surface currents also remain resolution
sensitive. A small energy residual is not an accurate formation result.

Half timestep changes pulse velocity below 0.0000710%; moving the magnetic
matching surface changes it below 0.000333%. Moving actual reflecting gas ends
from ±1 to ±2 m at a common matching extent of 3 m changes central |z|<0.5 m
velocity by 0.8174% / 0.8372%, and whole-domain source work by about 2.1%.
Spatial and material-boundary accuracy dominate those time/matching errors.
The change of gas ends is a physical-domain sensitivity, not an open-end model.

The nr=384 original shock with the strip completes to 0.02 s: B L1/B1 is
0.03065827, pressure L1 0.1753145 Pa, frozen-ratio defect 5.11e-13 and response
residual **1.01130e-7 J**, versus 0.00173278 J for the standard compatible case.
There are 985 steps and no additional adaptive limiting, but 1534 forced and
3064 neighboring affected cells. Minimum density is 0.0496321 kg/m³ without
a floor; fluid energy-flux residual is 3.28e-15 J. This preserves interior shock
accuracy and removes the observed growth of the strip's energy error through
nr=384. A separate nr=384 time study is not claimed.

All **32 new recorded cases** (12 pulses, 15 skin cases, four end controls and
the fine shock) are complete. No numerical jobs remain. Full fields, inputs,
source/environment fingerprints and the report are in
`outputs/boundary_verification` and `docs/BOUNDARY_VERIFICATION.md`.
An independent annular polynomial checks field restriction to 4.44e-16;
all saved-field hashes and the common frozen engine fingerprint match.
The final six-panel scientific figure was visually inspected after fixing
overlapping resolution-axis labels.

The next numerical implementation must address boundary/flow accuracy and
free-interface evolution using these source/skin/current checks. Keep the strip
experimental until its driven-field accuracy is adequate. Real gas ends,
transport/ionization/Hall validity, stronger shock/cold-flow and curved-equilibrium
accuracy, converged formation and a dynamic experimental UI remain in scope.
The Stage 2 formation gate and full objective remain unchanged.

## Moving material / infinite vacuum checkpoint

That checkpoint had **155 passing tests (61.39 s)**. That continuation was
verified progress toward the active objective, not completion. All 26 evolution
cases and 27 static evaluations completed. No numerical jobs remain live.
Source/configuration/environment and all lossless field hashes are checked;
both scientific figures were visually inspected, with logarithmic tick-label
collisions corrected. The full report is MOVING_PLASMA_VERIFICATION.md and the
derivation is MOVING_PLASMA.md, with records in outputs/moving_verification.

The new `MovingPlasma` reference has actual two-dimensional material motion,
adiabatic gas only inside its evolving surface, regular a=psi/r² induction
labels, analytic magnetic/pressure virtual work, an exactly curved matching
sphere and decaying magnetic modes extending to infinity. No exterior fluid,
pressure/density floors, numerical dissipation or energy repair are introduced.
This prototype is ideal-only with prescribed uniform bias; fixed-chamber
resistive references and production defaults remain unchanged.

Independent uniform-sphere/dipole energy and homologous force errors converge
approximately at second order. At nr=32 and scale=0.9, relative energy/force
errors are 0.000724094 / 0.002616710. Translating the sphere by 0.3 m exposes
spurious net axial force -0.00796587→-0.000122393 N at nr=4→32; it refines
but precludes a claim of exact magnetic axial-momentum conservation.

For exact free gas expansion at 0.1 s, velocity relative mass L2 errors are
0.0726360, 0.0268221, 0.00924936 and 0.00302204 at nr=4/8/16/32. Finest
surface-position error is 0.000138224 m and pressure-cell-average relative L2
error 0.000101608. Actual material mass is unchanged. The exact sphere radius
is 1.009950494 m; polygonal geometry and cell entropy projection are not repaired.

For the initialized spherical FRC, initial acceleration mass RMS decreases
0.713087→0.243325→0.0828970→0.0285246 m/s². Maximum artificial speed after
0.04 s decreases 0.0764241→0.0327842→0.0120209→0.00635480 m/s; finest surface
displacement is 0.000127825 m. This is measured equilibrium drift, not formation
or stability. Moving the matching sphere 2→3 m at nr=16 changes the artificial
velocity by at most 0.000331791 m/s RMS (12.2083% relative to the small reference
drift). Quadrature 5→7 changes it below 1.65e-12 m/s RMS; modes 4→8 below
2.15e-6 m/s. A tiny vacuum residual does not erase matching/mesh error.

Velocity Verlet and source-work integration show second-order time refinement
when the explicit cap controls actual steps. In the Gaussian uniform-bias ramp,
nr=8 response errors at dt=0.001/0.0005/0.00025 s are 0.000448658,
0.000112158 and 0.0000280390 J. However driven velocity changes **47.1666%**
between nr=8 and nr=16, and peak speed rises 0.166450→0.290958→0.524906 m/s
at nr=4/8/16. The driven ideal current-layer response is explicitly unconverged.
Source-work accuracy is not a reason to interpret it as physical formation.

Vertex psi is preserved near roundoff, while interpolated interior material
flux is not exact. A prescribed nonuniform deformation has relative flux errors
0.00386486→0.000929701→0.000229928→0.0000573222 at nr=4/8/16/32. Initialized
FRC interior flux drift at nr=32 is 1.85284e-6; source nr=16 is 0.000631995.
These are numerical interpolation effects, not physical resistivity or
reconnection. The regular-axis stretch approximation and low-energy/distorted
mesh modes remain explicit. No shock capture or remapping is provided yet.

The next resistive checkpoint below implements the scalar diffusion and heat
extension against the independent spherical-decay eigenmode. This ideal report
retains its historical source, including ideal source-flow limitations.

## Resistive moving material / vacuum checkpoint

The full suite has **170 passing tests (62.11 s)**. All 13 fixed-geometry decay
cases, independent stiff-mode controls and 20 actual moving cases completed.
Full configurations, histories, snapshots, source hashes and inspected figures
are in `outputs/resistive_moving_verification`; the method and report are
`docs/RESISTIVE_MOVING_NUMERICS.md` and
`docs/RESISTIVE_MOVING_VERIFICATION.md`. Frozen production/measurement SHA256 is
`452e51b4405b8a44fc6191f6716fa4e19e5798ba367bb84a0953ce3c6703503b`.
Report, analysis and plotting sources have separate recorded hashes.

`MaterialResistivity` uses a consistent plasma-only weighted mass and the same
full-space magnetic stiffness to solve scalar physical induction and recover
Jtheta. Implicit midpoint gives positive per-cell eta J² heat and a separate
source-work/energy check. `ResistiveMovingPlasma` atomically composes two
resistive half-steps with ideal material motion, updating both flux labels and
cell entropy. Fixed mass is retained and no exterior fluid, artificial
resistivity, pressure floor or energy repair is introduced. Eta=0 follows the
unchanged ideal evolution. Failure of the second resistive stage after trial
heating/motion leaves all accepted arrays, time and accounting unchanged.

An independent continuum conducting-sphere eigenmode matches a vacuum dipole
with continuous field at the surface. Maxwell/Faraday/Ohm derivatives and
all-space quadrature independently verify its current, energy and heat. At
t=0.1 s and nr=4/8/16/32/64, magnetic-field relative volume L2 errors decrease
0.0926223→0.00523351, current errors 0.0521209→0.000206685, and local-heat
errors 0.124056→0.000458581. Q1-gradient field error is about first order;
flux, current and cell heat are about second order. Axis-field and one-sided
surface-current errors are reported separately. A second radial eigenmode,
time, tolerance, matching-radius, mode-count and quadrature controls are complete.
The geometry is held fixed for this diffusion fixture; it is not force balanced.

A raw midpoint step of an independently diagonalized decaying spatial mode
produces a false amplitude -1/3 instead of exp(-4), despite a 2.22e-16 J energy
residual. Adaptive step doubling controls field, current and local heat, accepts
the unmodified two-half solution, and retains positive heat. Tolerances 1e-3
through 1e-6 reduce absolute modal-amplitude error 3.30e-4→4.13e-6 and preserve
the tested positive decay. This is local error control, not a global error,
positivity or arbitrarily small-field topology guarantee.

The actual moving cases start from the analytic FRC and then undergo resistive
relaxation, with or without a Gaussian uniform-bias ramp. Reversed topology is
initialized, not formed. With D=eta/mu0=0.1 m²/s, the driven nr=16→32 comparison
differs by 63.821% in mass-weighted velocity, 29.357% in magnetic change,
3.67623% in total current and 5.4924% in local heat density. Density change
differs by 54.0889%. These material-coordinate differences use coincident
material sampling and actual volume restriction, independently checked with
polynomial and geometric integration tests. They are not absolute solution
errors. Higher physical D=1 still gives 53.915% velocity difference; changing D
is a physical change, not numerical stabilization. The large flow error remains.

At nr=8, halving outer time caps gives approximately fourfold reduction in
velocity differences and response residuals. Driven energy residuals at
dt=0.001/0.0005/0.00025 s are 4.39132e-4/1.09776e-4/2.74435e-5 J; the last
velocity comparison differs by 0.0149697%. At nr=16, moving the magnetic
matching radius 2→3 m changes driven velocity by 0.539718%; mode/quadrature
changes are smaller. Tightening diffusion tolerance 1e-4→1e-5 leaves identical
accepted stages in these D=0.1 cases and therefore identical results. D=1 cases
do reject/subdivide diffusion trials. Mass is unchanged; maximum sampled
resistive-stage energy residual is 1.679e-14 J and vacuum relative residual
4.451e-13. These are distinct from full moving-state accuracy.

Next inspect spatial modes and surface thermal/force accuracy and add circular
coils at moving coordinates with the geometric external-source derivative.
An independent exact dipole/loop force fixture is specified in MOVING_COIL_PLAN.md.
Shock capture/remapping, material distortion, physical transport/ionization/Hall
regime, actual ends, converged formation and the experimental UI remain required.
All overarching issues and the Stage 2 entry remain open.

## Surface-resolution and thermal-compatibility checkpoint

This continuation implements optional `radial_grading` in `MovingDomain` and
records actual angular resolution in the resistive case contract. The default
uniform implementation is preserved: all archived arrays in the new source
nr=8/16/32 baselines are bitwise identical to the preceding checkpoint. The
frozen production/shared-measurement fingerprint is
`2749bf18bff1bdab93b1c338b7944ec7a6773d23c39076700f83a4faec9b90c8`.
All 21 dynamic and 13 independent cases completed, including the supplemental
finest quarter-time control (798 accepted motion steps). Analysis/render and
supplemental reproduction helpers are separately hashed. Full configurations,
lossless snapshots, histories, comparisons and figures are retained in
`outputs/surface_verification`; see SURFACE_NUMERICS and SURFACE_VERIFICATION.

Nine new independent tests bring the suite to **179 passing tests (80.30 s)**.
Revolved-polygon moments check graded geometry; exact free-gas expansion and
conducting-sphere decay check the actual kernels. Current/local heat remain
approximately second order under refinement. Grade-2 nr=32 current and local
heat relative errors are 0.000829269 and 0.00182061; Q1-gradient B error is
0.0105070. Grading trades interior resolution for surface resolution and does
not automatically improve the smooth free-gas reference.

Material-field sampling now inverts the actual initial sector geometry instead
of comparing equal normalized coordinates on different polygons. Exact convex
cell intersections resolve discontinuous thermal fields before positive Gaussian
integration. Independent affine field and geometric/indicator integrals verify
this diagnostic. Across all 18 comparisons, increasing integration order 4→6
changes reported relative norms by at most 3.38e-9 of the norm. This removes
the observed diagnostic-quadrature uncertainty without smoothing evolved fields.
The velocity norm integrates actual Q1 fields with initial physical mass;
thermal fields are compared pointwise rather than conservatively averaged onto
coarse cells. These definitions differ from the preceding report and are explicit.

The new velocity differences are 34.1273% for uniform nr=16/32 and 55.2917%
for grade-2 nr=16/32. Fine grade-2 angular doubling changes velocity 15.8362%.
Successive fine time controls change it 5.39874%, then 2.75800%; the ratio near
1.96 and nonmonotonic current differences do not establish asymptotic order.
Boundary-speed differences remain 0.04457785 then 0.01810780 m/s. The finest
quarter-time run has peak speed 0.5668651 m/s, minimum density 0.0555440 kg/m³,
physical heat 0.2738125 J and response residual 1.97602e-6 J. No floor, energy
repair or change of physical D is used. Spatial and time uncertainty remain.

Saved-state force decomposition exposes cancellation between large pressure
and magnetic surface accelerations. Recovered weak-current J cross B is a
separate discrete diagnostic, not an exact continuum reference. Mass remains
unchanged; maximum sampled resistive-stage energy and vacuum residuals are
1.386e-13 J and 7.566e-13. These do not establish local boundary accuracy.

The initial thermal-compatibility audit independently recovers the boundary
current from one-sided Maxwell derivatives with relative error below 4.7e-11.
The pressure equation gives positive equatorial initial p-dot of 3.75 Pa/s for
the unforced FRC and 1.379547904 Pa/s for the source case; a persistent smooth
p=0 trace under continuous-B/no-extra-stress assumptions would require zero.
The exact conditions and primary sources are in SURFACE_COMPATIBILITY. Keep
these initialized cases as stress tests, investigate a justified edge closure
or verified nonsmooth expansion, and do not substitute arbitrary eta tapers,
removed edge heat or easier fixtures for physical formation validation.

Both final scientific figures were visually inspected; colliding logarithmic
minor x-axis labels were corrected. The strict report serialization was fixed
to write diagnostic arrays as numeric lists, then the complete assembly was
rerun using all 18 verified comparison caches. All jobs are terminal. The goal
and relevant beads issues remain open with physical surface closure and the
independent spatial-coil force reference as concrete next work.

## Electron-current inertia and conditional transport checkpoint

This continuation adds `MaterialCurrent` with physical me and frozen prescribed
ne/eta profiles, counting current kinetic energy alongside magnetic energy and
positive per-cell collisional heating. Its sparse mixed Ampere/Ohm solve keeps
vacuum elliptic with no exterior carriers. Increment formulation preserves
stationary current; diagonal algebraic scaling improves axis-conditioned
roundoff without changing physics. A final exact-ideal-limit regression also
preserves tiny material flux under a much larger source-field change, with its
Ampere residual explicitly checked. Independent massless-limit and variable-
coefficient physical-volume energy checks pass. This is a separate linear
component, not a hidden change to the moving plasma or a complete edge model.

Thirteen new tests bring the full suite to **192 passed (79.10 s)**. Primary NRL
collision coefficients in SI agree with independent published cgs numbers.
Coulomb logs are never floored; invalid classical closure is explicitly
unavailable. Reference Stage 0/1 metadata now carries conditional thermal,
Hall/electron-inertia/Debye and drift diagnostics, retaining the original
prescribed conductivity and distinguishing inferred Drude rates. No absence of
a warning establishes validity; neutrals/ionization and tensors remain open.

All 35 current cases completed, including the added thin-layer nr=64 control.
Uniform inertial sphere decay, collisionless modified-Bessel screening and
an exact variable-density/resistivity profile independently test continuum
fields/current/heat. For de/R=.2, nr=8/16/32/64 current errors are
0.0249752/0.00768733/0.00247993/0.000828866; finest B and local-heat errors are
0.00523256 and 0.00104859. Finite inertia retains the weak derivative error in
the interpolated initial field, yielding roughly order 1.5 current convergence;
surface grading alone worsens this in the recorded control. This is explicit.

Strong screening de/R=.05 exposes false reversal on nr=8: minimum Bz/bias
-0.00330320 and 40.2668% of sampled material volume negative, while the analytic
field is positive everywhere inside. The amplitude is small relative to bias
but would invalidate a sign-only topology claim. Finer recorded samples have no
negative Bz; nr=64 current error is 0.00480722 and magnetic-field error 0.0730778.
No clipping, carrier floor, anomalous resistance or heat cutoff removes the
coarse artifact. Sample positivity is not a between-sample topology guarantee.

An independently diagonalized fixed spatial mode gives temporal halving ratios
4.002221 and 4.000555. Matching-radius, angular, quadrature, mode-count and
grading controls are retained. True-sphere energy quadrature orders 48/96 agree
within 9.172e-14 relative; maximum whole-run energy residual is 8.414e-15 of
the reported continuum energy scale and maximum equation residual 2.341e-12.
These do not erase current/field error or validate nonlinear formation.

The conditional scan records 43 assumed states and 41 prescribed-density-edge
samples. The original Ar+ case gives di/R=4.09721, de/R=0.0151831 and electron
mean-free-path/R=0.110877. Conditional perpendicular momentum and Z=1 parallel
resistivities are 8.57541e-5 and 4.37521e-5 ohm m; their limits are distinct and
neither is automatically installed. Fixed charge/temperature/ionization are
inputs, not equilibrium predictions. The prescribed vanishing-density edge
and its fixed current are a diagnostic illustration, not a thermal solution.

All source/configuration/lossless archive identities and companion documents
are checked. Frozen component measurement SHA256 is
`c66fb837ef6af44768c87d3a88e624f9ec75e318fd24a4cfd62a82752725908f`;
report, plotting, time audit and transport scan are separately hashed. Both
final figures are visually inspected. Full records are in
`outputs/current_verification`, equations in CURRENT_INERTIA and
TRANSPORT_REGIME, and results in CURRENT_VERIFICATION. All jobs are terminal.
EXTENDED_PLASMA_PLAN defines the next Hall/toroidal wave reference and the
remaining current-energy, thermal and edge integration requirements. The full
goal remains active; a component milestone is not FRC formation validation.

## Linear Hall / toroidal reference checkpoint

The full suite now has **208 passing tests (78.87 s)**, including 16 independent
Hall/finite-mass tests. `hall_parameters`, `hall_modes` and `hall_cylinder` add a
bounded uniform-density incompressible reference with regular poloidal/toroidal
fields and velocities. Exact center-of-mass/reduced-mass coefficients follow
parent species momentum; drift energy, bulk/magnetic exchange and positive heat
are paired. There is no exterior gas and no claim that its restricted periodic
reference boundary is a real plasma/vacuum or perfect-conductor boundary.

All 49 cylinder cases and 16 parallel modes (three time steps each) completed.
Actual cylinder field errors decrease from approximately 5.2% at nr=8 to 0.626%
at nr=64. Both helicities and prescribed scalar resistive damping are measured;
physical finite mass recovers ion/electron cyclotron asymptotes in the parallel
limit. A poloidal perturbation generates toroidal field and flow through Hall;
the MHD control does not. Physical volume integrals independently verify all
three energies, weak current and local midpoint heat; angular momentum is
conserved in this reference. The original argon target has not been replaced
by the separately labeled proton verification fixtures.

Parallel time-halving ratios are 3.998876–3.999721. Initial interpolation in the
slow cylindrical branch excites small fast secondary modes, giving nonmonotonic
apparent projection damping at coarse time steps. An independent full-spatial
helicity exponential and additional 2048/4096/8192-step controls resolve this:
relative time errors decrease 1.58885e-7 → 3.97212e-8 → 9.93030e-9. Remaining
spatial/projection error is explicitly separate. A one-step collisionless mode
still has 21.46% phase error despite roundoff energy conservation; no energy-only
accuracy claim is made.

Maximum whole-run energy-plus-heat residual is 6.63520e-15 relative, equation
backward error 2.38842e-17 and independent physical/modal energy disagreement
5.12417e-14 of initial energy. Cylinder source fingerprint is
`868ece2de1730a51d4aaf363652e3a19139c8289e046eda3188b57a0a45b0e3f`.
Full source/config/archive/history checks pass; records and inspected scientific
figures are in `outputs/hall_verification`. Equations and findings are in
LINEAR_HALL_NUMERICS and LINEAR_HALL_VERIFICATION. All jobs are terminal.

Next implement TOROIDAL_INTEGRATION_PLAN's conservative angular momentum and
toroidal induction in the bounded compressible core, verify geometric rotation/
pinch forces and nonlinear transport, then integrate Hall work and finite-mass
thermal/current transport from the parent species equations. The current
linear reference must not be inserted as an incomplete nonlinear closure.
Heated edge physics, spatial coils, actual vacuum/ends, shocks/remapping,
kinetic ordering, converged formation and experimental UI remain required.
The full goal remains active, and Stage 2 remains guarded.

## Nonlinear compressible toroidal checkpoint

The full suite passes **224 tests (85.76 s)**. Six new production modules add a
separate compressible state/solver with toroidal field, rotation, conservative
axial angular momentum and full magnetic/centrifugal stresses. Exact annular
r² moments enter pressure recovery and rotation/field energy, including the
axis. Both first/second-order zero-toroidal RHS recover the existing poloidal
solver, whose API and historical results remain unchanged. Failed final stages
leave accepted state, time, heat and limiter accounting untouched.

Scalar toroidal resistivity uses a positive compatible tridiagonal current mass
and matched Faraday/Poynting/heat work. It exactly recovers the independently
known 8Db diffusion of g=a+b r² at the axis and wall. Axis Ez is solved rather
than incorrectly set to zero; the conducting outer wall still has Ez=0. No
energy repair, current/density floor or hidden physical resistivity exists.
Massless MHD and fully filled periodic/conducting reference boundaries are
explicit; this is not nonlinear Hall or the free-vacuum formation interface.

All 33 declared cases completed. Exact rotation, pinch and combined equilibria
have maximum spurious radial RMS speed 3.99221e-16 m/s. A finite-amplitude
rotating entropy wave tests actual nonlinear angular-momentum/energy transport;
its nz=192 reconstructed density/rotation relative volume errors are
0.000126545/0.000125193. Torsional-wave reconstructed rotation error is
0.000462097; amplitude controls independently separate quadratic feedback from
spatial error. Stronger contrast and first-order diffusion controls remain
explicit, and no single stronger-contrast run is declared converged.

Independent conducting-Bessel decay reaches nr=128 (6689 accepted steps):
reconstructed toroidal field error 3.35910e-6, first-annulus g error magnitude
1.15482e-6 and physical heat relative error 2.07328e-5. Cell-average error is
separately recorded (1.95079e-6), not confused with field reconstruction error.
Physical quadrature orders 6/10 change reported relative norms by at most
8.21825e-17. Density time-halving ratios are 4.1723/3.9905 and angular momentum
3.7153/4.3801 against a finer same-grid reference, with limiter effects explicit.
Final mass/total-energy relative changes remain within 4.52346e-15/1.13086e-15.

All cases were regenerated under the final fixture/measurement contract, which
honors nondefault gamma and checks continuum fields between cell centers.
Frozen production/measurement SHA256 is
`24ce9b8cdb121167d78ab2ee75792bb9ace362695dd852cbc2384065f4dc69f2`.
Complete source/config/archive/history checks and scientific visual inspection
are recorded in TOROIDAL_MHD_VERIFICATION and `outputs/toroidal_verification`.
All solver, measurement, test and report jobs are terminal at this checkpoint.

Next implement NONLINEAR_HALL_PLAN: compatible nonlinear Hall induction and
Poynting work on the new compressible state, using independently derived
finite-sound-speed cylindrical modes and nonlinear transfer checks. Electron
pressure/thermal and finite-mass current transport, heated-edge compatibility,
spatial coil work, actual vacuum/chamber ends, shock/remapping robustness,
kinetic ordering, converged formation and the experiment UI remain required.
The full objective remains active and the original argon target is unchanged.

## Historical completion audit at the 224-test checkpoint

| Requirement | Authoritative evidence | Status |
|---|---|---|
| Conservative dynamic core | mhd and toroidal modules, independent tests, MHD/positivity and TOROIDAL_MHD_VERIFICATION reports | Bounded nonlinear poloidal/toroidal core implemented; positive traces without cell floors; full shock/update robustness and cold-flow accuracy remain bounded by tests |
| Physical source/plasma/vacuum coupling | Fixed-chamber, boundary, moving/resistive and SURFACE_VERIFICATION reports; SURFACE_COMPATIBILITY audit; Stage 2 guard still raises | Partial: moving resistive free surface and explicit surface refinement are implemented; conditional heated-edge incompatibility, unresolved flow/current, spatial coils, mesh robustness and actual ends remain |
| Regime-appropriate transport / Hall physics | Conditional Coulomb/local regime module, finite-current-inertia component and linear cylindrical Hall/toroidal reference; LINEAR_HALL_VERIFICATION and EXTENDED_PLASMA_PLAN | Partial: linear Hall/finite-mass waves and nonlinear compressible toroidal MHD verified; nonlinear Hall coupling, chosen formation closure, moving electron/ion thermal dynamics and physical edge/kinetic treatment remain incomplete |
| Full independent numerical verification | Prior dynamic/shock/vacuum/surface/current studies; conditional Coulomb checks; 49-case cylindrical Hall/16-mode parallel and 33-case nonlinear toroidal studies with independent species, transport, time and physical-field/energy verification | Partial; 224 passing tests; known false reversal on unresolved current layer, substantial thin-layer B error, moving-flow/cold-pressure/shock robustness and physical surface closure remain |
| Converged FRC formation | Only older stationary-conductor results, unconverged topology duration | Not demonstrated |
| Reproducible experimental interface | Stage 0/1 CLI/playback; dynamic and coupled benchmark scripts, histories and metrics | Partial; no coupled dynamic experiment UI |

Beads issues `frc-goal`, `frc-vacuum`, `frc-transport`, `frc-accuracy`, and
`frc-formation` preserve the remaining scope; `frc-mhd` records the completed
bounded reference core milestone. All numerical jobs at this checkpoint are complete. The optional request for target plasma parameters is not a blocker
to numerical verification or development.

The full objective remains active until all completion evidence is present.

## Nonlinear Hall spatial checkpoint: progress, failed axis gate retained

The full suite passes **244 tests (84.28 s)**. New independently checked
finite-sound-speed Hall modes derive all three branches from compressible
momentum, pressure and Faraday; original electron/ion and Maxwell residuals,
physical energy, conducting/axis traces and limiting dispersions are verified.
The declared reduction has massless cold electrons and adiabatic ion pressure.
A separate finite-amplitude variable-density analytic snapshot verifies local
Hall current, induction and Poynting identities before judging numerical work.

The experimental HallCurl/HallInduction/HallPoynting components cancel assembled
Hall magnetic work through a common electric quadrature and an exact adjoint
of the existing magnetic energy. Both poloidal and toroidal rates are present;
there is no carrier floor, Hall Joule heat, extra fluid force or energy repair.
They are **not enabled in ToroidalMHD**, whose prior verified path is unchanged.

All 42 instantaneous spatial cases completed. Nonlinear nr=8/16/32/64/128
physical magnetic-rate errors decrease 0.259461/0.132568/0.0666923/0.0333967/
0.0167042, and current error reaches 0.00235369. However, first-annulus gdot
error approaches 0.214581 instead of zero. A shrinking physical axis volume
hides that coefficient defect in global norms. This explicitly fails the Hall
integration gate; global conservation is not accepted as physical accuracy.
Local energy-transfer defect decreases to 0.000109583 relative but is not
installed as heating or used to repair magnetic energy.

Maximum assembled Hall work ratio is 2.72060e-17. Independent archive energy
replay backward error is 2.60690e-17 using the actual summation scale. Reference
quadrature 10/14 changes cell comparisons by at most 4.75452e-16; common
quadrature 5/7 changes the nr32 physical induction error by 1.97065e-15.
All three finite-pressure branch/sign Hall-rate errors refine to 3.52–4.00%
at nr32. These are spatial rates, not time-evolved waves. Direct quadratic
Hall amplitude powers are 1.999999987 and 2.000000001; the ill-conditioned
opposite-full-field subtraction and its small-amplitude failure remain saved.

Frozen production/measurement SHA256:
`35cea1122c521983ac6808a8d8f6d27c3cf430428561e2f0de45819598a93996`.
Every declared case/config/source/archive identity passed audit; both final
scientific figures were visually inspected. All test, measurement and report
jobs are terminal. Details: COMPRESSIBLE_HALL_REFERENCE, HALL_WORK_NUMERICS,
HALL_WORK_VERIFICATION and `outputs/hall_work_verification`.

The next action is to fix regular-axis polynomial current/electric/energy
compatibility, then local Poynting transfer and stable nonlinear Hall time
evolution. The original SSPRK2 fluid method cannot simply receive an undamped
imaginary-spectrum Hall RHS. Full electron/ion thermal and finite-mass transport,
heated-edge/vacuum/real-end closure, spatial coils, shock/remapping robustness,
kinetic ordering, converged formation and interactive experiments remain open.
No formation claim or original-argon target change is made. The full goal is
active; this is progress rather than completion or an external blocker.

## Compatible nonlinear Hall checkpoint: verified progress

The full suite passes **257 tests (88.52 s)**. The corrected electric/magnetic
spaces commute under physical Faraday in x=r²/2,z, with explicit compatible
magnetic energy carried through pressure recovery, fluid stages and diagnostics.
The old cell-energy path and failed Hall operator remain reproducible controls.
There is no fitted axis coefficient, density/current floor, energy reset or
extra Hall force. A nonlinear implicit midpoint couples Hall and positive
resistivity to ideal compressible fluid half steps; failed trials are atomic.

All **54 main plus 9 supplemental cases** completed. Axis gdot error falls
0.0203290→0.00526712→0.00132551→0.000331905→0.0000830089 for nr8→128, replacing
the historical nonzero 0.214581 limit. Uniform/strong-density, zero-component,
quadrature and radial-only controls retain convergence. Global Hall work ratio
is at most 4.41337e-18; independent physical energy replay is within 7.94556e-16
of its absolute work scale. Local energy-transfer defect reaches 9.59915e-5
relative and converges near second order, but is not zero or repaired as heat.

The actual coupled solver reproduces all three finite-pressure Hall branches.
At nr32, physical magnetic errors are 3.9960%/4.4048%/3.9631%; initial field
representation already contributes 3.9807%/4.3985%/3.9631%. This is a remaining
first-order spatial field limitation, not accurate-everywhere B inferred from
small budgets. Slow/intermediate/fast velocity errors are
1.0118%/0.2128%/0.2759%. Slow duration is 0.5 radian, the others a full period.
Slow mode retained amplitude improves 0.6961→0.9679→0.9969, so coarse damping
is substantial and explicit. A 100-fold smaller amplitude hardly changes its
nr16 error. Both propagation signs, Hall-free, nonlinear amplitude and solver
tolerance controls are archived; physical and coefficient errors stay separate.

Fixed-grid time-halving magnetic error ratios are 3.94368/3.98028 for the wave
and 3.94105/4.02325 for nonlinear evolution. The finite-amplitude 1-microsecond
flow comparison against restricted nr64 has nr32 magnetic/momentum/density
differences 0.2432%/0.3791%/1.2872% of the fine reference's **change**, rather
than its large background. This does not make nr64 an exact continuum oracle.
Independent conducting Bessel decay has nr64 field/axis/physical-heat errors
1.21612e-5/9.87096e-6/1.33088e-4, with time controls separating spatial error.

All evolved histories conserve total energy/mass to at most
6.63912e-16/5.40031e-16 relative. Independent initial/final energy replay, actual
terminal times, positivity, physical heat monotonicity, source/protocol identity
and lossless archive hashes pass audit. Angular momentum uses an absolute
normalization where its mean is zero. All three scientific figures were visually
inspected. Frozen engine/reference/measurement source is
`9c14e7c757b6d3a51bc0988faff5ff704cd58e2abb8b5ebaada646de7e4d4dd8`;
supplemental protocol and report hashes are separately recorded. All numerical,
report and test jobs are terminal. Details: COMPATIBLE_HALL_NUMERICS,
COMPATIBLE_HALL_VERIFICATION and both `outputs/compatible_hall_*` directories.

### Current completion audit and next action

| Requirement | Current status |
|---|---|
| Nonlinear compressible Hall core | Implemented and independently verified in a bounded conducting periodic reference; first-order physical B representation and limited shock/update scope remain explicit |
| Electron/ion thermal and finite-mass closure | Incomplete: cold massless electrons, single represented gas, prescribed eta; physical heating partition, electron pressure/enthalpy and nonlinear drift-energy terms still required |
| Plasma-vacuum, edge, ends and sources | Earlier moving/free-vacuum references retained; heated-edge incompatibility, spatial coils, real ends and unresolved flow/mesh errors remain |
| Converged FRC formation | Not demonstrated; initialized equilibria and waves are verification fixtures |
| Interactive physical experiments | Existing Stage0/1 interface only; guarded Stage2 and dynamic experiment UI remain required |

Next implement TWO_TEMPERATURE_PLAN's **discrete conservative charge flux and
electron pressure/enthalpy pairing**, then warm electron/ion evolution and
independent collision exchange. The weak work current is not automatically a
pointwise divergence-free cell particle flux; its dual identity must be mapped
and verified before thermal advection. The electron enthalpy drift is required
alongside pressure Ohm/Poynting terms. Finite-mass convective current, drift
energy and momentum stresses must follow from the same species equations.

Return thereafter to original heated-edge/vacuum/real-end physics, spatial coil
work, shock/remapping robustness, kinetic ordering and converged formation.
No easier or colder fixture replaces the original argon target. Beads goal,
transport, vacuum, accuracy and formation issues stay open. The full objective
is **active**; this checkpoint is progress, not completion or an external block.
