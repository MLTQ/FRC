# Smooth finite-mass transport: improvement and remaining magnetic failure

The separate smooth reference improves shared mass/enthalpy and electron entropy
transport. The full suite passes **336 tests in 118.11 s**. All **31 rate cases
and 27 trajectories** completed and passed the independent audit. These are bounded verification fixtures;
accurate formation, a heated plasma–vacuum surface and physical transport remain
unfinished. Conservation and positive species states do not prove local accuracy.

The equations, boundaries, electric-kernel representation and absence of added
diffusion are specified in [SMOOTH_FINITE_MASS_NUMERICS](SMOOTH_FINITE_MASS_NUMERICS.md).
The scientific figure was visually inspected after correcting overlapping
logarithmic tick labels. The previous limited method and its source/data are retained unchanged in
[FINITE_MASS_VERIFICATION](FINITE_MASS_VERIFICATION.md).

![Pressure improvement and unresolved magnetic response](../outputs/smooth_finite_mass_verification/accuracy.png)

## Protocol and independent references

The study repeats all 21 earlier rates and 18 trajectories with identical physical
inputs, grids and time windows. It adds two hundredfold weaker wave-rate controls,
eight magnetized translating-contact rates (four grids, smooth/limited entropy),
three full-transit unmagnetized contacts, a limited-entropy comparison, reverse
flow, half dt, and three short magnetized contacts. Full case definitions,
parameters, initial/final fields and every accepted diagnostic remain in the
records. Previous slow/medium wave branches still have rate evidence only.

The exact contact is independently specified by n(z−Vt), pe(z−Vt), pi=p0−pe,
uniform bulk velocity and uniform axial B. Its full curl-free Ez=−pe_z/(ne)
balances the separate species pressure forces. Tests independently differentiate
both parent momentum/adiabatic equations and the full energy flux. There is no
forcing or heating in these contacts. Physical argon/electron masses, n0=1e20
m^-3, R=.035 m, L=.24 m, V=2400 m/s, gamma_e=5/3 and gamma_i=1.4 are used.
Density and pe modulation amplitudes are .25 and .12; reference temperatures
are 5 and 1 eV. Unmagnetized cases run 100 microseconds, one whole transit.
The B0=.005 T cases run only 50 ns at a maximum dt of 1 ns; they must not be
presented as a full-transit magnetized demonstration.

Initial nodal pe is exact. Initial cell pi complements the actual Q1 pe average
to maintain discrete constant total pressure. That creates measured initial
pi/energy projection error; evolved states are never reset to the reference.
The nonlinear manufactured parent-species and exact-mass linear wave oracles
are the independent references described in the preceding verification document.

## Complete-rate accuracy

Errors are relative L2 against physical rates, not equilibrium backgrounds.
Nodal electron pressure uses dual-volume weights; ion pressure compares cell
means. Full physical B_t/J_t/E norms are integrated independently.

| nr | Manufactured mass rate | Manufactured pe rate | Manufactured pi rate | Slow-wave pi rate |
|---:|---:|---:|---:|---:|
| 8 | .000605757 | .556797 | 5.54737 | 4.41005 |
| 16 | .0000382054 | .183055 | 1.43950 | .610186 |
| 32 | .00000239326 | .0614256 | .371728 | .0691434 |
| 64 | .000000149664 | .0210460 | .0979971 | .00883534 |
| 128 | .00000000935531 | .00731328 | .0268603 | not run |

The finest manufactured pi maximum decreases from 7426.95 to 3225.44 Pa/s
at nr64→128, unlike the old increasing maximum. At nr128, its first-annulus
relative error is .0104550. Radial/axial/angular momentum rate errors are
2.98846e-5/1.35953e-7/7.07279e-6; total-energy rate error is .000175310.
The original n0=1e20 fixture improves less: nr64 pi rate error is .0857483.
The coarse nr8 manufactured pi rate is slightly worse than the prior method;
the improvement is not universal across every grid and quantity.

Slow nr32 pi rate improves from 117.248 times its signal to .0691434;
fast nr32 pi rate improves from .454781 to .0213591. Weak-amplitude controls
retain those thermal improvements. The fast axial-momentum rate error changes
.0573952→.00590473→.00151359 under amplitude factors 1,.1,.01: nonlinear
finite-amplitude truncation contributes to that particular linear-oracle error.

Contact pe-rate errors at nz16/32/64/128 are .00395518/.000282350/
.0000182693/.00000115188. Their limited-entropy counterparts are
.293983/.0410243/.0278047/.0186435. Ion-rate errors refine approximately at
order two because their initial species representations differ.

## Evolved pressure and contact accuracy

Fast-wave errors are normalized by their physical oscillation, including
continuous magnetic and species-pressure representation error.

| nr | Magnetic field | Electron pressure | Ion pressure | Prior ion pressure |
|---:|---:|---:|---:|---:|
| 8 | .158332 | .129221 | .472967 | 15.4730 |
| 16 | .0795212 | .0340228 | .147025 | 3.24947 |
| 32 | .0398052 | .00876142 | .0538315 | .595516 |

The nr32 weak-wave pi error is .0538235. Magnetic accuracy is unchanged by
the shared mechanical transport correction. In the nonlinear manufactured
nr32 run, ion error falls only from .00362760 to .00325299 of full pressure,
and from .783580 to **.702663 of its actual change**. The original-density
counterpart is .00244043/.527145. Small errors against a large pressure
background therefore still conceal inadequate resolution of physical changes.

After a full unmagnetized contact transit:

| nz | Nodal pe error | Cell-mean pi error | Total pressure error | Spurious velocity RMS (m/s) |
|---:|---:|---:|---:|---:|
| 16 | .00149711 | .00821875 | 1.28395e-5 | .0561255 |
| 32 | .000103524 | .00136420 | 4.54902e-7 | .00115827 |
| 64 | .00000663411 | .000316986 | 3.82396e-8 | .0000256064 |

The limited-entropy nz32 contact gives pe/pi errors .00839935/.0369084
and spurious speed .919900 m/s. Reversed flow gives pe/pi errors
.000103739/.00130932. Half dt gives .000103524/.00136420; timestep
reduction is not removing the spatial/projection error. These are distinct
measured solutions, not forced reflection or equality checks.

These nodal/cell-mean values are **not continuous physical field errors**.
Separate `contact_physical` audits integrate Q1 electron and cellwise ion
pressure/density, and their pressure/(e n) temperatures, over physical volume.
Initial continuous pe/pi errors at nz64 are .0000743749/.0110716; initial
electron/ion temperature errors are .000496586/.0127338. Spatial representation
alone therefore creates appreciable ion-field error before any time integration.

After the full transit, continuous pe/pi errors at nz16/32/64 are
.00188155/.0449188, .000319681/.0221656, and .0000750472/.0110717.
The finest electron/ion temperature errors are .000498286/.0127339.
These are relative full-field norms; the report also retains dimensional and
sampled maximum errors. The B0=0 contact's false magnetic RMS is
8.34e-12/2.95e-12/1.14e-11 T, with current RMS
.00182757/.000782152/.00405749 A/m². This small numerical drift is **not
monotonically convergent** and remains explicit. The much shorter magnetized
contacts have final B−B0 RMS below 1.51e-18 T and J RMS below 1.03e-9 A/m²;
their short windows are not comparable long-time stability evidence.

An independent closed-form check verifies these initial pressure quadratures.
Let h=pi/nz, s=sin(h)/h, and A be the pe modulation amplitude. Mean-square pe
error is A²[(2+cos(2h))/6+1/2−s²]; mean-square pi error is
A²[cos²(h)/2+1/2−cos(h)s]. Relative norms divide by the corresponding exact
mean-square species pressure. The formulas agree with physical quadrature at
nz16/32/64; results are preserved in `contact_norm_control.json`.

Because a full transit returns the exact contact to its initial state, the exact
pressure change is zero. The report retains dimensional error in that change
and a null relative error, rather than dividing by numerical roundoff. Magnetized
contacts separately retain nonzero small-change and absolute B−B0/J drift norms.

## Conservation, time and entropy audit

The largest accepted-history energy-budget residual is 3.36e-16 relative to
initial total energy. Normalized conserved-integral drifts are at most 4.97e-15;
near-zero net momentum is normalized to noncancelled physical input scales.
All accepted densities and species pressures are positive. Every source archive,
configuration, field hash, saved input, history endpoint and reported norm is
checked. Physical component energy/particles/entropy are independently integrated;
parent total-energy change checks manufactured forcing. Quadrature refinement
changes reported relative L2 norms by at most 3.09e-10 absolute.
Independent physical energy replay differs from stored energy by at most
2.55e-14 relative; parent source-work discrepancy is at most 5.08e-12 J.
Minimum accepted pe/pi/n are .734719 Pa/.137941 Pa/9.18491e17 m^-3.

Manufactured temporal difference ratios are 15.978/15.992 for B,
16.200/16.141 for pe and 15.593/15.835 for pi. These agree with fourth-order
time convergence in that smooth fixture. Fast-wave B/J/pe ratios are
15.478/15.475/15.261, while pi is .7957 with successive dimensional differences
1.67e-14 and 2.09e-14 Pa. A universal fourth-order pressure claim is therefore
not supported. The tiny time differences do not explain the much larger
spatial errors. Full finite-mass energy subtraction limits very weak pi signals.

Unforced nonlinear physical electron log-entropy change per particle is
−5.77e-6/−1.63e-6/−4.18e-7 at nr8/16/32. Ion values are
.000253247/.0000489541/.0000113619. Contact nz32 forward/reverse electron
changes are +1.19e-6/−1.22e-6; ion changes are +4.88e-5/−4.83e-5.
Exact smooth adiabatic contacts conserve these quantities. These measured
defects disprove any inference of a discrete entropy inequality from conserved
total energy; physical shock partition remains unverified.

## Unresolved slow magnetic response

The pressure correction is not a successful slow-wave solution. Its physical
magnetic-rate relative errors are **1004.98, 263.069, 66.8838 and 16.8431**
at nr8/16/32/64. Corresponding current-rate errors are 1194.34, 302.992,
74.2153 and 18.1646. Slow conserved-energy-rate errors remain 106.207/26.5782
at nr32/64. Exact magnetic and total-energy-rate norms are comparable; this
cannot be dismissed as a tiny energy denominator alone. A hundredfold
weaker amplitude retains the nr32 magnetic error.

The bulk-motion and Hall electric terms each have about 19200 times their net
magnitude. Errors in the lower-order nodal velocity and weak-current fields
spoil this cancellation. A separate linear control uses the original physical
density/masses/bias, with perturbation-only input to exclude stationary-bias
roundoff in the diagnosis. It substitutes analytic fields only in the forcing:

| nr | Both numerical | Exact velocity | Exact current | Both exact |
|---:|---:|---:|---:|---:|
| 8 | 1004.98 | 204.151 | 1184.71 | .158539 |
| 16 | 263.069 | 51.2841 | 309.221 | .0795489 |
| 32 | 66.8838 | 12.8374 | 78.5564 | .0398095 |
| 64 | 16.8431 | 3.21043 | 19.7755 | .0199092 |

These are diagnosis controls, not adopted corrections or alternative physical
initial conditions. Even the both-exact control retains first-order physical
field representation error. A consistent higher-order velocity/current and
electromagnetic representation is the next required derivation and verification.

## Electric kernel and retained failures

The original contact rejection is preserved with its exact input, forcing and
source identity. A large curl-free electric field has coefficient norm 61.4,
versus a non-kernel remainder of 2.20e-10. The new algebraic separation retains
the full E for work and strict relative residuals for the small-signal solve.
Independent original mixed equations verify weak signals beside a large kernel.

Original full-equation backward errors are below 9e-17 in the preserved contact
audit. However, recombined full-E Faraday replay still has **.016679 relative
error against the tiny rate**. This is explicitly retained, not renamed as a
small relative residual. The exact contact magnetic rate is zero; the measured
roundoff-sized rate is not evidence of physical induction accuracy. Separable
kernel extraction is not a general decomposition and may not cure other
ill-conditioned fields.

The first smooth study completed 31 rates and 18 trajectories before a NumPy
boolean caused contact JSON serialization to fail. That attempt and two orphan
contact field archives remain in `smooth_finite_mass_initial_attempt`. Only the
boolean serialization was corrected; all cases were regenerated with a new
source identity. Missing original histories are never described as verified.
The audit confirms exact reproduction of all 49 old complete records and their
arrays, and both orphan contact field archives, excluding only elapsed time and
the documented source identity. All fresh contact records have full histories.

Known driver limitation for future work: actual step failures preserve the last
accepted state and history, but the outer RuntimeError/ValueError handler is too
broad and could mislabel an unexpected post-integration metric failure as an
initialization failure. No such failure occurred in these 27 verified records.
Narrow that handler before expanding the measurement protocol; retain partial
states and distinguish solver failures from reporting failures.

## Reproducibility and next work

Final measurement source:
`dd060015dc510d312f22a53588d2809ad806de02de1dadb4db956e8887ab4ead`.
Initial smooth attempt:
`5c045bea6b32e7467d99b9366c373c51ca2ba2e73f3ef8439d2830222ed5bd97`.
Both include complete source snapshots. The old finite-mass trajectory source
`3edae52928f5cdaa20ca240f11d6fa729a137be327c9b1ffc64f09c6440f42fb`
and earlier rate/moving/current identities remain unchanged.

Fresh data: `python -m examples.smooth_finite_mass_measurements --kind rates`
or `--kind trajectories`, with an explicit new `--directory`. Full read-only
audit: `python -m examples.smooth_finite_mass_report`. Exact configurations,
field hashes, source/import closures, accepted times, budgets, analytic norms,
independent physical energy/entropy and report helper snapshots are retained.
`checkpoint.json` records final report/figure/source/test/document hashes and
terminal process status. No remote repository or publication is configured.

Next correct the shared velocity/current magnetic cancellation while retaining
the improved transport, independent pressure/energy tests and actual physical
masses. Then continue collisions/heating partition, conduction/ionization,
the original heated free edge with exact vacuum, real ends/coils, shocks and
kinetic validity, converged formation and the experimental interface. No 3-D
stability or particle-confinement conclusion follows from these fixtures.
