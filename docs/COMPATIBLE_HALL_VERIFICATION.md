# Compatible nonlinear Hall verification

**257 tests pass (88.52 s); all 63 declared measurement cases and the report audit
complete.** The previous 21.46% axis defect is corrected by commuting physical
electric/magnetic spaces and an explicit compatible magnetic energy. A new
compressible solver now evolves nonlinear Hall/resistive induction, density,
flow and total energy together. This is progress toward the full goal, not a
formation result or a physically complete transport/edge model.

The method, reduction and boundaries are in COMPATIBLE_HALL_NUMERICS. The
main 54 cases contain 23 instantaneous spatial audits and 31 actual coupled
evolutions. Nine separate controls add two slow-wave amplitude tests and seven
resistive Bessel decays. Inputs, initial/final arrays, every accepted diagnostic
time, actual step counts, versions, individual records and collection records
are archived in `outputs/compatible_hall_verification` and
`outputs/compatible_hall_controls`. No numerical job remains live.

## Regular axis, independent nonlinear snapshot and work

The unchanged analytic HallManufactured reference has all magnetic/current
components, nonuniform positive carriers and finite Hall nonlinearity. Its
corner psi and cell g are integrated independently at orders 10/14. The new
operator uses electric quadrature 5 and fixed magnetic mass quadrature 16.

| Radial cells (nz=2nr) | First-annulus gdot relative error | Physical induction L2 error | Current L2 error | Local energy-transfer defect |
|---:|---:|---:|---:|---:|
| 8 | 0.0203290 | 0.258067 | 0.0790832 | 0.0228900 |
| 16 | 0.00526712 | 0.132229 | 0.0240889 | 0.00604306 |
| 32 | 0.00132551 | 0.0666313 | 0.00993017 | 0.00152994 |
| 64 | 0.000331905 | 0.0333869 | 0.00473307 | 0.000383670 |
| 128 | 0.0000830089 | 0.0167028 | 0.00234160 | 0.0000959915 |

The axis coefficient converges near second order rather than being hidden by a
shrinking cell volume. Uniform carriers, contrast 0.8, zero toroidal field,
zero poloidal perturbation, and radial-only refinement at fixed nz=256 retain
convergence. At nr32 the axis errors with quadratures 3/5/7 are
0.00132754/0.00132551/0.00132551. Independent reference quadrature changes are
at floating-point scale. There is no axis multiplier, carrier floor or repaired
field/energy. Physical field-rate/current norms remain partly first order;
the piecewise field spaces are not made higher order by coefficient convergence.

Maximum assembled Hall work ratio is **4.41337e-18**, normalized by actual
absolute work. Independent 24-point physical B dot Bdot integration replays
global cancellation to **7.94556e-16** of that work scale and local magnetic
rate to **8.58683e-16** relative. Poynting and magnetic local rates do not agree
exactly; their measured defect above is retained, not converted into physical
heat or used to reset state energy. The older failed data are verified against
their original commit f28714c before comparison, never relabeled as current.

![Axis and spatial refinement](../outputs/compatible_hall_verification/spatial.png)

## Actual compressible three-branch propagation

The independent finite-pressure continuum eigenproblem uses R=0.035 m,
L=0.24 m, ne=1e20 m^-3, Ar ion mass, B0=0.005 T, p=16 Pa and gamma=5/3.
The explicit reduction has cold massless electrons and adiabatic represented
gas. Frequencies are 2794.26185, 218192.58861 and 760901.20611 rad/s.
The slow branch advances **0.5 radians**, not a full cycle; the other branches
advance 2pi. The baseline velocity eigenvector normalization is 1 m/s.
Both axial propagation signs, fast amplitudes 10/100 m/s, and a Hall-free
compressible mode are included. All histories use the actual coupled solver.

| Branch | Radial cells | Physical B perturbation error | Physical velocity error | g cell error |
|---|---:|---:|---:|---:|
| Slow | 8 | 0.325930 | 0.334357 | 0.279706 |
| Slow | 16 | 0.0868162 | 0.0630078 | 0.0320121 |
| Slow | 32 | 0.0399602 | 0.0101183 | 0.00426753 |
| Intermediate | 8 | 0.188891 | 0.0609239 | 0.0873846 |
| Intermediate | 16 | 0.0886046 | 0.00963196 | 0.0124623 |
| Intermediate | 32 | 0.0440484 | 0.00212781 | 0.00257264 |
| Fast | 8 | 0.158172 | 0.0494019 | 0.00843516 |
| Fast | 16 | 0.0792156 | 0.0109568 | 0.00143723 |
| Fast | 32 | 0.0396311 | 0.00275904 | 0.000335784 |

Physical norms integrate reconstructed fields and perturbations between cell
centers, excluding uniform bias from the denominator. Initial physical magnetic
errors at nr32 are already 0.0398069/0.0439855/0.0396306. These mostly explain
the final magnetic error floor and show why small g cell error is not the same
as accurate physical B everywhere. Initial pressure error is also retained:
the true continuum total energy is averaged without fitting pressure to the
discrete magnetic/kinetic convention. For the fast nr32 wave, initial/final
pressure perturbation errors are 0.00569522/0.00489492.

An independent two-template toroidal-mode projection separates damping and
phase. For the slow branch, retained fundamental amplitudes at nr8/16/32 are
0.696088/0.967888/0.996909, with phase errors
-0.0154754/-0.00132704/+0.0000962448 rad. The coarse result is substantially
damped. For intermediate/fast nr32, amplitudes are 0.997434/0.999889 and phase
errors 0.00144522/0.000310638 rad. Residual outside the fitted fundamental is
recorded; projection does not erase it from field errors.

Slow-wave nr16 amplitude reduction 1→0.1→0.01 m/s changes magnetic error
0.0868162→0.0868143→0.0868142 and velocity error
0.0630078→0.0630074→0.0630073. Thus the observed coarse slow-wave error is
predominantly numerical, not nonlinear physical feedback. Opposite axial
propagation errors agree to the recorded numerical precision. Finite-amplitude
fast controls remain close to the linear reference but are not exact linear
solutions or independent finite-amplitude oracles.

![Evolved waves and initial representation error](../outputs/compatible_hall_verification/waves.png)

## Time, nonlinear redistribution and resistance

Fixed-grid time studies compare full archived fields with a finer **same-grid**
trajectory. Norms divide by that reference's change from its own initial state,
so uniform bias/pressure cannot conceal the time error. Actual step counts
match all prescribed counts. On nr8, fast-wave phase 1 with 8/16/32 steps
against 256 steps gives magnetic errors
0.00125633/0.000318568/0.0000800368: halving ratios **3.94368/3.98028**.
Nonlinear end time 3e-7 s with 4/8/16 against 128 steps gives
0.00176887/0.000448832/0.000111559: ratios **3.94105/4.02325**.
Momentum ratios are 4.03266/3.81894 for the wave and 3.95359/4.03645 for
the nonlinear case. Density/energy controls are independently recorded.
An additional small-mesh test compares midpoint directly against an independent
DOP853 integration of the same semidiscrete nonlinear ODE.

The finite-amplitude HallManufactured initial state is unbalanced and evolves
for 1e-6 s with full fluid feedback. Against the nr64 conservative/corner state
restricted by exact annular volumes, nr8/16/32 relative differences are:

| Variable | nr8 | nr16 | nr32 |
|---|---:|---:|---:|
| Magnetic norm | 0.0429466 | 0.0110370 | 0.00243185 |
| Momentum norm | 0.0952850 | 0.0244391 | 0.00379106 |
| Density norm | 0.163483 | 0.0501868 | 0.0128724 |
| Total energy norm | 0.0465119 | 0.0136882 | 0.00258008 |

These denominators are the restricted fine-reference **change**, not its total
state. This is a nested conservative-state comparison, not proof that nr64 is
an exact continuum trajectory. Peak speed increases 176.05→179.68→180.91→181.55
m/s; nr64 minimum pressure is 15.2732 Pa. The stronger-density and eta=1e-4
controls complete, but single controlled runs are not themselves convergence
studies. The resistive nonlinear run diagnoses 2.41740e-4 J physical heat.

Conducting Bessel decay uses Btheta=C J1(qr) exp(-eta q² t/mu0), qR a J0 zero
so Jz(R)=0, C=1e-7 T, eta=1e-4 ohm m, exponent 0.2 and 64 steps. The exact
continuum magnetic loss is 3.26536e-13 J. Coupled fluid feedback is O(C²), not
identically zero. Magnetic/heat errors converge quadratically:

| Radial cells | g weighted L2 error | First-annulus relative error | Physical heat relative error |
|---:|---:|---:|---:|
| 8 | 0.000759511 | 0.000574449 | 0.00843721 |
| 16 | 0.000192536 | 0.000156451 | 0.00211760 |
| 32 | 0.0000483337 | 0.0000397655 | 0.000530303 |
| 64 | 0.0000121612 | 0.00000987096 | 0.000133088 |

At fixed nr16, 32/64/128/256 steps give heat errors
0.00211960/0.00211760/0.00211711/0.00211698: spatial error dominates there.
Pressure change is not falsely equated with Ohmic heat, since compressional
work is also O(C²). No incompatible resistive J1-wall wave is used.

![Time, nonlinear mesh and decay controls](../outputs/compatible_hall_verification/controls.png)

## Conservation, provenance and remaining limits

Across all evolved histories, maximum relative total-energy drift is
**6.63912e-16**, mass drift **5.40031e-16**, and absolute axial angular-momentum
drift **2.58494e-26 kg m²/s**. Zero-mean angular momentum is not used as a
relative denominator. Decay histories have fewer columns, so their mass/angular
checks use archived endpoints. Initial/final mass, kinetic/magnetic/internal
energy and angular momentum are independently replayed from physical fields;
every accepted step reaches the declared time with positive pressure/density
and monotonic nonnegative physical heat. Hall-only heat is exactly zero.
Conducting/periodic/axis traces and finite residuals are audited. Final face
divergence is at most 7.5634e-16 T/m in the main evolution set.

The nr16 fast run with tolerances 1e-8 and 1e-10 takes identical iterations;
versus 1e-12 its magnetic trajectory difference is 3.80129e-7 of the reference
change. The stricter reference is not declared exact. Solver tolerance, fluid
time stability, Hall phase accuracy and spatial accuracy remain separate.

Frozen engine/reference/measurement SHA256:
`9c14e7c757b6d3a51bc0988faff5ff704cd58e2abb8b5ebaada646de7e4d4dd8`.
The supplemental control script has its own hash in every configuration. Report
sources and figure hashes are separate in `report.json`. The audit validates
exact protocol coverage, unchanged source, individual/collection JSON agreement,
NPZ hashes, finite arrays, mesh coordinates, all history lengths and terminal
times. All three figures were visually inspected. Reproduce with:

```sh
.venv/bin/python -m examples.compatible_hall_measurements
.venv/bin/python -m examples.compatible_hall_controls
.venv/bin/python -m examples.verify_compatible_hall
```

The solver is not enabled as the coil-driven Stage 2 formation model. Separate
electron/ion thermal and finite-mass transport, physically justified heated
edge/vacuum/end conditions, spatial coil force/work, shock/remapping robustness,
kinetic ordering, converged formation and a dynamic experimental UI remain
required. Next: TWO_TEMPERATURE_PLAN's discrete charge/enthalpy/pressure work
pairing, then independent warm-species evolution. The original argon target
and all earlier unresolved moving-interface/flow errors remain in scope.
