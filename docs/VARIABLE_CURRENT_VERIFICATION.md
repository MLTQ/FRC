# Changing-density compatible-field verification

This study verifies an **instantaneous variable-density induction/energy
component** and measures an explicit lumped pressure-work bridge. It does not
evolve the full moving plasma, replace the earlier smooth solver, or demonstrate
FRC formation. The [numerical contract](VARIABLE_CURRENT_NUMERICS.md) states the
density representation, energy derivative, source and boundary assumptions.

The final protocol has **14 independently audited cases** and **351 passing
tests in 123.17 s**. Source is
`fc106dc19b9bc208aaf7d9b36a43dc69ea49f4d49d126b4f4dfd2ecd85c59e97`, with
raw arrays, source/companion archive, audit, report and figure in
`outputs/variable_current_verification`. All workers are terminal.

The initial complete protocol and audit are preserved in
`outputs/variable_current_before_range_guard` under source
`310809e22d2b4f1aa3f0e11fbe28ec72bf19f724fc57295eb4338e3fb999fb01`.
Its 351-test suite passed in121.09 s. Final review added explicit rejection of
reconstruction overflow and unrepresentable physical inertia, without a density
floor. **Every saved numerical array in all14 cases reproduces bit for bit**
after that guard; `range_guard_replay.json` records the comparison.

![Measured component convergence](../outputs/variable_current_verification/verification.png)

## Physical component accuracy

The default manufactured cylinder uses physical Ar+/electron masses, n0=1e18
m−3, axial bias 0.005 T, and smooth radial/axial density and fields. Density rate
comes from actual integrated prescribed particle flux through faces. No fitted
current, effective mass, frozen n, numerical atmosphere or energy repair is used.
An analytic source forces only the declared instantaneous induction component;
this is not claimed as a complete parent-species trajectory.

Relative physical L2 errors for degree 5, with twice as many axial cells:

| Radial cells | B_t | J_t | Local current-energy rate | Density contribution to that rate |
|---:|---:|---:|---:|---:|
| 8 | 1.58106e-5 | 1.47994e-5 | 2.73321e-5 | 3.60406e-5 |
| 16 | 5.05338e-7 | 4.72892e-7 | 8.99917e-7 | 1.18018e-6 |
| 32 | 1.58811e-8 | 1.48605e-8 | 2.84902e-8 | 3.73408e-8 |
| 64 | 4.96985e-10 | 4.65041e-10 | 8.93213e-10 | 1.17057e-9 |

Finest B_t axis-relative error is 2.15568e-11, RMS 7.09560e-9 T/s and sampled
maximum 2.37572e-8 T/s. J_t axis-relative error is 5.44158e-10, RMS 0.642381
A m−2 s−1 and sampled maximum 3.75616 A m−2 s−1. These absolute values and the
large manufactured rate's relative errors are retained separately.

Degree-3 B_t errors are .0147796/.00199058/.000243538/.0000293840. Additional
nr32 cases vary quadrature8/10/12, n0 from1e16 to1e20 m−3, density contrast, and
the sign of the prescribed continuity flow. Degree-5 B_t remains about1.58811e-8
in these controls. The sparse iteration count grows to1419 at n0=1e16; the
physical electron mass and strict1e-10 solver tolerance are unchanged.

## Thermal transport remains the limiting component

The new Hessian-derived graph flux pairs exactly with the new Ampere current.
The previous diagonal-g stream does not: its electron relative net-work mismatch
is .0325780/.00798787/.00198728/.000496217 over the same grids. The corrected
pressure bridge retains near-roundoff lumped thermal work and fixed-n logarithmic
entropy, **but still uses a low-order thermal representation**.

| Radial cells | Electron pressure E, physical L2 | Ion pressure E, physical L2 | Electron U_t, lumped nodal L2 | Ion U_t, lumped nodal L2 |
|---:|---:|---:|---:|---:|
| 8 | .00821914 | .00722727 | .343915 | .344153 |
| 16 | .00205379 | .00179833 | .222161 | .220097 |
| 32 | .000514248 | .000450095 | .123670 | .122369 |
| 64 | .000128628 | .000112574 | .0635129 | .0628264 |

The pressure electric field converges roughly at order2. The lumped nodal
entropy-rate error is about6.35%/6.28% on the finest grid, with much slower
convergence than the field component. These nodal norms compare nodal rates to
continuum point values with the declared lumped volume weights; they are not
misrepresented as continuous high-order thermal quadrature. Both species are
checked separately before multiplication by their physical mass fractions.

The radial error profile exposes a stronger local limitation. The first four
radial-node rows contain **99.8072% of the electron squared error** at nr64,
peaking at row1. Their local relative L2 error is .310637/.373928/.387672/.390882
over nr8/16/32/64: it does **not** converge in that shrinking near-axis region.
At nr64 the interior relative error is .000976709 and the last-four wall region
is .0332307. The whole-domain improvement partly reflects the shrinking volume
of the unresolved axis region. The new `thermal_error_regions` report records
the full radial profile and absolute RMS/max; this is not hidden by a global
energy or convergence claim. The pressure bridge is unsuitable for adoption as
an accurate high-order thermal update in its current form.

## Independent work and archive checks

Eight focused tests use analytic polynomial fields and density, exact annular
moments, physical energy integrals, central directional derivatives, an
independent original three-block solve, and complex-step continuum curl/Faraday/
particle-continuity checks. Invalid positive-cell reconstructions reject without
modifying inputs. The final guard additionally rejects finite coefficients whose
reconstruction overflows and density too small for finite positive inertia.

Independent root-product field and cumulative-primitive density reconstructions
at quadrature12/14/16 audit the saved arrays. Maximum measured defects:

| Quantity | Maximum |
|---|---:|
| Local energy/rate agreement with independent quadrature14 | 4.56198e-15 |
| Density-energy pullback vs independent physical work | 2.16217e-15 |
| Independent quadrature12 vs16 local integrals | 1.01147e-15 |
| Full original-equation backward error | 9.93891e-13 |
| Solver's checked residual | 1.32456e-11 |
| Full-array Faraday error relative to rate | 5.76072e-10 |
| Magnetic/current/density energy exchange backward error | 8.19274e-16 |
| Surface-flux continuity vs analytic volume rate | 3.65973e-15 |
| Integral continuity backward error | 3.62756e-18 |
| Species pressure-work backward error | 1.14084e-15 |
| Fixed-n logarithmic entropy backward error | 9.21608e-17 |

The full-array rate-relative Faraday error is distinct from the physical
energy-norm solve and original-equation backward errors; it is explicitly not
claimed below1e-10. Local harmonic-bias cross work is retained, even though its
global integral is zero. Small net rates do not supply the sole error scale.

## Next action

Derive a common higher-order thermal density/energy/transport representation
using the same physical density reconstruction and its adjoint. Resolve the
nonconverging near-axis lumped entropy rate before coupling pressure into full
moving evolution; mesh refinement alone is insufficient for the local error.
Then add the velocity/vorticity and bulk/current stress, enthalpy and Poynting
transport with the complete conservative energy partition. Repeat independent
parent-species manufactured rates and complete wave/contact trajectories with
small-pressure-change, axis/max, time/quadrature and atomic rejection checks.

The original heated edge/exact vacuum, transport/ionization, chamber/end/coil
coupling, shocks/entropy, kinetic validity, converged formation and interactive
experiments remain open. Prior smooth source `dd060015...` and compatible linear
source `8cfb1a2f...` remain unchanged.

Reproduce in a fresh output directory with `python -m
examples.variable_current_verification <directory>`, then run
`examples.variable_current_audit` and `examples.variable_current_report` with
the same directory. Use the repository virtual environment and a writable
Matplotlib cache. The focused check is
`python -m pytest tests/test_variable_compatible_current.py -q`.
