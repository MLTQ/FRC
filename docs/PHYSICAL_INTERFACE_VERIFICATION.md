# Local physical budgets and the larger-grid accuracy failure

This checkpoint continues the complete FRC sandbox goal. It is a filled,
periodic, fully ionized argon verification fixture, with physical electron and
ion masses. It is not formation, an experiment design, or a confinement result.
The previous matched-response checkpoint is preserved. No production evolution
source, physical mass, acceptance tolerance, or accepted state was changed.

## Original local equations

Read PHYSICAL_INTERFACE_AUDIT_DESIGN.md. An independent exact one-sided sampler
now integrates particle, current, both species entropy, and total-energy fluxes
on every physical cell face. It retains both traces, their signed and absolute
jumps, every cell's actual volume rate, its outward flux, and their residual.
The periodic seam is an internal face; the axis has exactly zero face area.
No small boundary displacement, trace averaging, fitted flux or energy repair
is used. Current has no evolved charge density: its actual divergence and
normal-current jumps are measured against the full physical current gradient.

The original total energy flux is the sum of each species' kinetic and
enthalpy flux and E cross B / mu0. All physical masses, electron/ion drift
velocities, current inertia, actual electric field and harmonic bias are kept.
Independent manufactured smooth/discontinuous trace tests, an expanded
mass-moment energy-flux identity, and an analytic local polynomial transport
balance verify the audit. The last fixture is explicitly an instantaneous
transport identity, not a complete finite-mass Ohm solution.

Eight actual full-RHS initial continuum perturbations were measured, with
quadrature16/20 and complete fresh bitwise replay. Here F denotes fluid
elements/axial points and B denotes field radial/axial cells. Values below are
relative cell-balance L2 defects; unlike-unit channels are not summed.

| F | B | particles | current | electron entropy | ion entropy | total energy |
|---|---|---:|---:|---:|---:|---:|
| 4/32 | 8/32 | 1.26793e-2 | 1.07224e-3 | .367589 | 4.17265e-3 | .604850 |
| 4/64 | 16/64 | 2.95080e-4 | 8.87575e-5 | .0249428 | 9.87285e-5 | .0453847 |
| 4/128 | 32/128 | 6.04965e-6 | 5.98929e-6 | 1.56772e-3 | 1.95663e-6 | 2.81076e-3 |
| 8/128 | 32/128 | 5.14120e-6 | 5.98929e-6 | 1.56754e-3 | 1.72854e-6 | 2.81026e-3 |
| 16/128 | 64/128 | 5.14012e-6 | 5.99536e-6 | 1.56834e-3 | 1.72826e-6 | 2.81167e-3 |
| 8/128 | 32/256 | 5.16609e-6 | 3.79094e-7 | 1.35710e-4 | 1.73776e-6 | 2.69221e-4 |
| 8/256 | 32/256 | 1.30603e-7 | 3.79094e-7 | 9.74110e-5 | 4.02227e-8 | 1.73539e-4 |

The eighth case repeats F16/128,B64/128 at epsilon1e-5 rather than1e-4;
energy defect is2.81140e-3. Refining only radial meshes at axial128 does not
remove the axial defect. Refining field axial cells alone to256 helps, and
refining fluid axial points separately helps further. The final energy defect
is about16.19 times smaller than F8/128,B32/128. This is a measured consistency
result, not closure of the full nonlinear interface/formation gate.

The largest quadrature16/20 relative discrepancy among these cell L2 norms is
9.64e-10. Independently replayed axis/wall-element norms are retained: final
energy defects1.74812e-4/1.72022e-4. Absolute energy trace-jump integrals fall
from5.18584e-4 W to1.84378e-7 W. The coarse net exterior/global power is about
1e-17 W despite its .605 local defect: cancellation cannot certify local accuracy.
The exact cell divergence theorem and geometric face-flux accounting are saved
separately. Sources and all raw cell/face arrays are in
outputs/matched_interface_budgets and outputs/matched_interface_axial256,
with separate region, replay and verification archives.

## Full response and time controls

Both new coarse120-step and larger-grid60-step signed trajectories cover the
full100 microseconds, at epsilon +/-1e-4. Every original accepted state/stage
is saved. All four full evolution reruns reproduce every accepted state/stage
and non-runtime metric bitwise. All physical audits, histories and comparisons
also replay bitwise. Every job is terminal;457 tests pass in254.80s.

At120 coarse steps, final relative errors n/Br/Bz/J/E are
.000970955/.00188784/.000875620/.00274345/.00848749.
The60/120 physical differences Br/J/E are
.000213545/.000202453/.00101175. Along with30/60 differences these **do not
demonstrate fourth-order temporal convergence**. This gate remains open.

At larger F4/64,B16/64, final relative errors are:

| channel | initial | final |
|---|---:|---:|
| density | 1.69e-8 | .000817825 |
| radial velocity | 4.30e-9 | .000151659 |
| radial B | 4.10972e-5 | .000426294 |
| axial B | 3.73560e-6 | .00496073 |
| current | 2.77310e-5 | **.0431611** |
| electric field | .000274455 | .00575621 |
| electron pressure | 1.71452e-8 | .000743703 |
| ion pressure | 1.69314e-8 | .000627370 |

The intended component projections improve to about2.40216–2.40218, close to
the independent2.402175837. Nevertheless current error grows from the coarse
.00274791 to .0431611. **This is a failed full-response mesh-accuracy gate.**
Good projection, positivity, equation residual or global energy cannot replace
the physical-field comparison. No formation adoption or mesh convergence is
claimed. Complete saved physical histories retain parallel and orthogonal
content at every accepted endpoint, with quadrature12/16 and an independent
endpoint comparison to the frozen quadrature16 audit.

## Fast response diagnosis

The larger-grid current's orthogonal norm divided by the initial reference
norm rises from .00136713 at40us to .10368056 at100us. Subtracting the known
slow exponential between40/80 and60/100us windows isolates a diagnostic fast
direction without changing any evolved state. Its current grows by4.1939704
over20us, equivalent to71682.39/s, with .00533 window-to-window shape error.

The actual stationary full spatial RHS, evaluated directly without a time
integrator at both signs of that normalized direction, projects to71426.05/s.
Its current rate has .04651 nonparallel relative residual. Changing the probe
amplitude1e-4 to1e-5 changes the full current rate by1.57e-6; independent
quadrature16/20 agrees. This is strong evidence for a spatial contribution,
but is not a fully converged eigenpair or proof of a unique cause.

About97.83% of extracted current squared norm is azimuthal and98.07% of B norm
is axial. Axial Fourier index1 dominates; only2.58e-10 of current norm lies
above index16. Radial Bessel indices15/16 dominate Bz/density on a field with
16 radial cells. The first16 radial Bessel modes capture only~63% of density
and71% of Bz;32 modes capture~94% of both. Uncaptured fractions remain explicit.

The independently reconstructed strong curl(B)/mu0 differs from weak current
by .599823 in relative L2. This strong volume comparison excludes tangential
interface distributions; it must not be mistaken for a complete weak Ampere
identity. Applying the independent original continuum operator to the same
projected physical direction yields large radial-momentum/induction rate
discrepancies. The direction is near the radial resolution limit, not a
validated physical instability. Archives: matched_mode_extracted_response and
matched_mode_extracted_continuum, each with complete bitwise replay.

Exact subdivision of the same fluid tangent to8/16 elements preserves physical
fields within7e-14 while leaving the field mesh fixed. Its instantaneous
current rate changes by less than8.4e-8 relative, but the velocity rate changes
by factors45.2/75.2 of the original rate norm and becomes almost perpendicular
to the original velocity direction. Thus fluid refinement alone is not an
established repair or converged eigenmode. This separate control is fully
replayed in matched_mode_extracted_fluid_refinement. No filter/damping was added.

A separate representative uniform spectrum exposes a further time-resolution
risk. At60 steps, maximum omega*dt is118.27 on the coarse grid and318.04 on the
larger one. Gauss stability does not imply phase accuracy. This is a uniform
operator, not a bound on the nonuniform spectrum or a proof that a particular
fast branch was excited; it does not explain away the direct spatial growth.

## Small-signal energy scale

For U(n,C)=A*n^gamma*exp(C/n), a=delta n/n and
b=(delta C-s*delta n)/n, the thermal half-Hessian is

    U2 = U/2 * [(b+(gamma-1)*a)^2 + (gamma-1)*a^2].

Together with physical bulk kinetic, current inertia and magnetic quadratic
energies, this gives an initial positive normalization about v=J=0. Independent
70-digit differentiation of the original conservative EOS and an analytic
fixed-C second derivative verify it. It is **not a separately conserved free
energy** on the nonuniform growing background.

Initial scale at epsilon1e-4 is1.6526365e-11 J. Actual continuous total-energy
drift divided by that scale is at most6.42e-7 over the five paired audits,
including the new120-step and larger-grid runs. This is a more informative
normalization than background energy, but does not certify local or field
accuracy. All signs, quadratures, components and physical fingerprints replay.

## Required continuation

Keep frc-7mm.7/.6/.4 and the full FRC goal open. Resolve the radial near-grid
branch using independently matched radial/field/fluid spaces and fixed-physical-
wavenumber consistency controls, then repeat complete signed trajectories with
genuine mesh, time, amplitude and residual refinement. Do not suppress the
branch or use a fitted production eigenmode as a success fixture. Keep local
volume/face/jump diagnostics and small-signal energy checks in those controls.

Then complete heated moving exact vacuum, transport/ionization validity,
physical ends/coils, converged prescribed-waveform formation, and the experiment
UI. The user's~1m quartz tube,10–20cm diameter,0.1–100mTorr neutral fill, selectable
gas and Helmholtz-like bias remain exploration inputs, not recommendations.
The six-panel PNG/PDF in outputs/physical_interface_report is visually reviewed.
