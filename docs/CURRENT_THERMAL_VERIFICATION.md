# Finite-mass current and thermal verification

This is a bounded, fixed-positive-density, constrained-zero-barycentric-flow
reference. Physical electron mass, nonlinear current inertia and both species
thermal transport are included. Moving plasma, a free edge and FRC formation
are not established. Equations and boundaries are in CURRENT_THERMAL_NUMERICS.
The final full test suite passes **316 tests in 103.41 seconds**.

## Audited baseline

The 41-case collection has **40 completed trajectories and one explicitly
rejected massless trajectory**. Its independent archive/history/physical replay
passes, including the failure record; this is not a claim that every case passes
an accuracy criterion. Both generated figures were visually inspected.

Baseline source fingerprint:
`6c3e50a2d595ee1e2bbe14c5a6679d8da0855ab99143fa802f3422543dec8cd8`.
The exact pre-fix source is preserved in
`outputs/current_thermal_verification/source_snapshot.tar.gz`, SHA256
`cc5f892e1ba27c63c8164de7656c0741d5571b9c1cf754a49b739f094ed81f88`.
The collection/report, individual configurations, accepted histories and 41
initial/final field archives retain source and content hashes. The snapshot
also preserves the earlier report helpers and tests. Later precision changes
do not relabel these measurements as new-source runs.

Independent physical quadrature replays magnetic energy, reduced-mass current
energy and both species thermal energies. Maximum all-history relative drift
is 5.80e-16 for their sum and 3.87e-16 for independently evolved conservative
cell energy. Maximum electron/ion log-entropy changes are 1.65e-14/6.87e-16 per
particle. Local agreement of the two energy accounts remains a measured error.

## Original nonbarotropic fixture

The reference retains n0=1e20 m^-3, R=0.035 m, L=0.24 m, Te0=5 eV, initial
ion pressure 16 Pa, its original toroidal field, zero poloidal bias and final
time 1 microsecond. No easier imposed bias or density replaced this fixture.
The physical finite-mass trajectories use 80 equal time steps.

| Radial cells | Local pressure-equivalent RMS defect, Pa | Maximum, Pa | Minimum inferred ion pressure, Pa |
|---:|---:|---:|---:|
| 8 | 0.0591249 | 0.422662 | 15.5774 |
| 16 | 0.0158876 | 0.116989 | 15.8830 |
| 32 | 0.00404562 | 0.0309969 | 15.9690 |
| 64 | 0.00101606 | 0.00784713 | 15.9922 |
| 128 | 0.000254304 | 0.00196778 | 15.9981 |

The finest first-annulus maximum is 0.00178157 Pa. Doubling the finest run to
160 steps gives RMS 0.000254306 Pa and maximum 0.00196780 Pa. This particular
defect is spatially dominated and approaches second order. Success does not
establish a general positivity theorem. The conservative ion partition is
checked independently of the positive evolved species entropies; neither
energy account is reset to match the other.

With consistently massless electrons the nr=128 joint run rejects its 60th
trial at 0.75 microseconds after accepting 0.7375 microseconds. Its nonlinear
residual is 1.69355e-9 after 80 iterations, exceeding 1e-10. At the last accepted
state the maximum local defect is 1.04326 Pa and minimum inferred ion pressure
is **positive**, 14.9567 Pa. This is a nonlinear rejection, distinct from the
older pressure-only split's negative-ion-pressure failures near 0.95 microseconds.

## Independent nonlinear waves and time refinement

Exact finite-amplitude helical Bessel waves satisfy both species momentum
equations, Faraday and the full local total-energy law, with an explicit
mass-proportional constraint reaction whose summed work vanishes. Both
helicities, reversed axial propagation, target and inertia-dominated densities,
strong amplitudes and a massless control are retained. Physical interpolation
is independently evaluated, with errors normalized to the perturbation.

For n=1e18 m^-3 and vector-potential amplitude 1e-6 T m, magnetic L2 errors at
nr=8/16/32 are 0.158454, 0.0795402 and 0.0398093. Current L2 errors are
0.0831394, 0.0409820 and 0.0203837; first-annulus current error at nr=32 is
0.0283100. These roughly first-order physical field errors remain significant,
even though global energy is conserved. Strong target-density waves give
similar convergence; their nr=32 radial/axial face-power errors are
0.00382125/0.00421943. Face power includes Poynting, both species enthalpies
and nonlinear drift kinetic transport.

At fixed mesh, wave magnetic RMS differences from a 320-step reference are
3.73492e-8, 9.32040e-9, 2.30334e-9 and 5.48452e-10 T at 10/20/40/80 steps.
Their ratios 4.007/4.046/4.200 support second-order time convergence; the
reference has finite time error. Original-fixture time, quadrature, density-unit
scale and solver-tolerance controls are independently retained in `report.json`.

![Baseline convergence](../outputs/current_thermal_verification/current_thermal_convergence.png)
![Baseline limits](../outputs/current_thermal_verification/current_thermal_limits.png)

## Continuum and nonlinear-inertia limits

The local massless pure-toroidal principal-symbol analysis is independent of
the discrete solver. On the original fixture, the joint fixed-density symbol
has minimum discriminant -7.39664e7 m^2/s^2, corresponding to frozen short-wave
growth 4300.19 |k| s^-1. Retaining moving ions and density in the massless
first-derivative symbol still produces imaginary characteristic speed up to
4916.37 m/s in the scanned directions. Complex-step differentiation of the
continuum equations checks the symbol. This is a local warning, not a theorem
about the global nonlinear evolution or proof of the cause of an axis failure.

The independent non-Beltrami inertia forcing audit deliberately tests terms
that cancel in the exact helical wave. Volume error of the intermediate B*
forcing decreases 0.126911 → 0.00771241 across nr=8→128, but first-annulus
relative error stalls near 0.235399. The maximum error at nr=128 is 15.1558
T/(m s), versus maximum reference 66.0026. This defect remains recorded in
`outputs/current_inertia_spatial`; it must not be confused with physical Bdot
after the finite-inertia inverse. A separate independent Fourier/Chebyshev
physical-rate study distinguishes those two quantities.

That physical-rate study is complete and audited in
`outputs/current_rate_schur_verification`. It solves the same static mixed
equations through a symmetric positive magnetic-space elimination, checking
the reconstructed Ohm/Ampere/Faraday equations. Fresh direct-solve comparisons
pass on nr=4/8/16; saved earlier direct results on nr=8/16/32/64 agree within
5.71e-13 relatively. The largest reconstructed residual is 1.87e-13. This
auxiliary solve changes neither production time integration nor physical mass.

Independent Fourier/Chebyshev references use three radial orders 64/128/256.
Radial reference differences are below 7e-14; separate axial-resolution and
solver-tolerance controls differ by at most 3.21e-12 in cell means.

| Radial cells | Full physical rate L2 | Full first annulus | Convective-only L2 | Convective-only first annulus |
|---:|---:|---:|---:|---:|
| 8 | 0.0253926 | 0.0244104 | 0.109806 | 0.217729 |
| 16 | 0.00639101 | 0.00570131 | 0.0374250 | 0.142305 |
| 32 | 0.00160233 | 0.00126512 | 0.00913541 | 0.0675810 |
| 64 | 0.000401131 | 0.000281948 | 0.00168899 | 0.0236863 |
| 128 | 0.000100341 | 0.0000656197 | 0.000274914 | 0.00689978 |

The physical inversion substantially reduces the intermediate forcing's axis
defect. It does not erase that intermediate error or validate moving nonlinear
evolution. The isolated convective first-annulus error remains about 0.690% at
nr=128. The new figure was visually inspected. The earlier direct study's local
archives and remote partial log remain under `outputs/current_rate_verification`;
the costly remote nr=128 factorization was explicitly stopped after the
equivalent study completed. It is not labeled as a completed direct run.

![Physical-rate verification](../outputs/current_rate_schur_verification/physical_current_rate.png)

## Stationary-background precision correction

The frozen baseline rejects all four uniform 0.005 T controls at nr=8/16 and
dt=1e-8/1e-6 s. Exact continuum current is zero; weak matrix currents are only
1.32e-9 and 1.31e-8 A/m^2. Recomputed full currents were subtracted to form the
inertia residual, magnifying roundoff relative to the tiny physical signal.
Failures at time zero and unchanged inputs are preserved in
`outputs/current_bias_failure`, rather than overwritten by successful reruns.

The correction applies the same linear current map directly to the Faraday
increment. It checks Ohm on that small-signal scale and separately checks the
stored state's Ampere/Faraday backward error. The physical equations, electron
mass, boundary conditions, nonlinear tolerance and failure policy are unchanged.
No current is thresholded away, and no diffusivity, atmosphere or energy repair
is added. New-source repeated stationary and tiny-wave evidence is recorded
separately in `outputs/current_precision_verification/report.json`.

The correction protocol completes and audits **18 stationary controls and 15
evolutions**, under source
`80fd45bddda51871f7f36682e8d0bdc4495d0c41430191c25ef03e06ed47c4d2`.
Stationary cases cover ten steps each on nr=8/16/32, Bz=0.005/-0.005/0.5 T and
dt=1e-8/1e-6 s. Their largest field error divided by bias is 3.54e-14, and both
species entropies stay bitwise unchanged. The 15 evolutions include six original
fixture/base-wave reruns, six weak waves, and strong-bias/reversed-helicity/
strong-inertia controls. Maximum component-energy drift is 3.09e-16; maximum
stored-state and nonlinear residuals are 7.88e-17 and 1.95e-11 respectively.

At amplitude 1e-12 T m, the physical magnetic L2 errors are
0.1585613 → 0.0795541 → 0.0398110 at nr=8/16/32; current errors are
0.0833629 → 0.0410087 → 0.0203869. Amplitude 1e-9 T m gives nearly identical
relative errors. Thus the accepted solution retains tiny waves instead of
silently discarding them. Across the six repeated baseline cases, maximum
change in reported relative field/current/face-power errors is 4.49e-14;
the maximum local pressure defect changes by at most 1.14e-13 Pa.

For tiny waves, recovering a rate by subtracting two stored background fields
can lose significant relative digits. The independent audit retains that
relative discrepancy and its absolute state-roundoff bound. Passing the bound
does not claim precision beyond the represented state; continuum wave errors
still use the perturbation itself as denominator.
The strongest-background control's recovered Faraday rate loses up to 2.22e-4
relatively while remaining within its explicit absolute storage bound. Its
pressure-replay uncertainty is bounded by 1.56e-8 Pa, exceeding its tiny
~8.25e-10 Pa reported defect. That defect should therefore be interpreted as
roundoff, not physical precision at the quoted digits. Both limitations are
preserved in the independent report.

## Reproduction and next work

For the historical baseline, extract its source snapshot to a separate working
directory, use the versions recorded in its provenance, and run
`examples.verify_current_thermal` against a copy of the frozen data there.
Do not overwrite the baseline with current-source measurements. Run
`examples.verify_current_precision` in the current checkout for the correction
protocol. Both numerical paths reject source changes during measurement.
For the complete physical-rate study use
`python -m examples.current_rate_measurements --method schur --directory outputs/current_rate_schur_verification`,
then call `examples.verify_current_rate.verify` with that directory. Preserve
the old partial direct-solve archives separately.

The next physical step is an independently verified moving finite-mass stage:
drift stress in momentum, variable-density current energy, bulk/current
advection, velocity-gradient terms and the full energy flux must be derived
together. A manufactured moving, varying-density reference should first check
both parent species equations and local conservation. Updating only an inertia
matrix in the existing split is insufficient. Collisions/heating partition,
conduction/ionization, exact vacuum and the heated free edge, real ends/coils,
shock entropy/partition, kinetic validity and converged formation remain open.
