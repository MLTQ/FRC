# Coupled moving-plasma checkpoint (2026-10-04)

The actual regular fluid, current transport, separate species entropy/pressure
and magnetic equations now run together in a bounded filled periodic cylinder.
This is verified progress toward the experimental sandbox, not an FRC formation
or laboratory prediction. Read REGULAR_PLASMA_NUMERICS.md for equations,
boundaries, physical energy and the absence of artificial diffusion/heat.

## Reproducible evidence

The final regression suite passes **425 tests in 219.19 seconds**. All new local
study, replay, energy, refinement, reporting and test jobs have terminated.
No remote job was launched, no remote is configured and no publication occurred.

The main study has **25 cases**, all arrays and non-time metrics replayed bitwise:
16 original parent-rate comparisons, seven full fast-wave periods and two full
separate stationary spectra with nonlinear controls. Independent reconstruction
uses physical quadrature orders 16 and 20. Source (70 Python files):
`19f68fb08f949c4526bc1347d18652fb1eeb22d720d483a2074d6a71f7f40358`.
Raw arrays, source archive and terminal audit are in
`outputs/regular_plasma_verification/`. Its first nineteen records explicitly
link the preceding source; every one was fully recomputed under the final
source during audit. No result was relabeled without that recalculation.

Additional frozen evidence:

| Study | Evidence | Source SHA-256 |
|---|---|---|
| Six slow-branch amplitude/mesh controls, bitwise plus quadrature 16/20 | regular_plasma_slow_controls | 543ed10530781378499d4b650a75c8246e68f89b5fd0619d2cc541279d1b98ac |
| Seven continuous-energy audits, both quadratures replayed bitwise | regular_plasma_energy | e98f0d45d077f30f31945d0851a63e93d573c88dafe20e07c84da672e5ebc45a |
| Captured original strict rejection, exact input replay | regular_plasma_kernel_failure | 96b0c95e74c35d4a0499596a2f054bc38f8571498d14eaef4cdb4e2c29ea138a |
| Refined full coordinate matrices and complete unchanged wave | regular_plasma_refinement | bc562813c4df25bb032e693b4e59adf0782aa9e92ee29f1233a36b34a1a2cd3f |

These names are directories under outputs/. The initial coordinate-difference
attempt/source `8fc86897cbaa07dfb72a896226a70bc8de44e7ae80951a11b575d2430fa173d8`
is retained, including its terminal failed case 19; its intended case 20 was
not run after that failure. Initial direct-background energy and pre-input
failure captures are retained separately. Seventeen earlier source closures
remain unchanged. The report is `outputs/regular_plasma_diagnostics/convergence`
as PNG/PDF; it has been visually inspected.

## Original parent rates

Physical argon/electron masses are unchanged. The nonlinear contact has variable
density, unequal complementary pressures and axial speed 2400 m/s, with nonzero
curl-free E and zero exact current. At 4/8/16 degree-five radial elements:

| Relative physical L2 error | 4 | 8 | 16 |
|---|---:|---:|---:|
| Density rate | 1.56571e-6 | 2.55830e-8 | 4.04333e-10 |
| Electron pressure rate | 2.17819e-6 | 3.57247e-8 | 5.65129e-10 |
| Ion pressure rate | 3.09363e-5 | 5.59953e-7 | 9.10035e-9 |
| Full electric field | 1.04658e-2 | 4.99936e-4 | 1.82313e-5 |

Contact acceleration remains at roundoff relative to velocity/transit-time.
Residual reconstructed B/J rates are recorded as absolute errors, not silently
set to zero. This is an instantaneous comparison, not a completed contact transit.

For the fast compressible finite-mass wave, finest (16 element) relative rates
are velocity 1.66544e-6, E 1.41611e-6, B 2.31526e-6 and J 2.20055e-6. Density
and pressure rates are around 2e-9. Axis, wall and maximum norms are retained
in the raw records. The largest original induction residual over all rate
cases is 1.51037e-11 at the unchanged 1e-10 tolerance.

## Full unforced fast-wave evolution

The continuum period is 1.16517746482e-7 s. All actual fluid/entropy and magnetic
variables advance together, with physical current inertia and no source,
temperature reset, diffusion or energy correction. At 120 steps per period:

| Final perturbation-relative L2 error | 4 elements | 8 elements | 16 elements |
|---|---:|---:|---:|
| Density / both pressures (approximately) | 2.87864e-2 | 5.54700e-3 | 5.6764e-4 |
| Velocity | 1.05388e-3 | 1.07331e-4 | 1.52793e-5 |
| Magnetic field | 2.22398e-3 | 7.31870e-5 | 2.35089e-6 |
| Current | 2.15809e-3 | 7.41538e-5 | 2.59599e-6 |

The finest first/last-element density errors are 7.18346e-5/1.41226e-3; velocity
2.16760e-5/9.94362e-5; magnetic 5.35887e-7/6.94799e-7; current
3.19748e-6/5.47023e-6. Initial representation errors are separately saved.
No norm uses the large uniform density/pressure/bias as the wave denominator.

For 60/120/240/480 steps, successive final magnetic differences are 5.90301e-6,
3.68951e-7, 2.30596e-8; ratios **15.999437 / 15.999863** verify fourth-order
time convergence. Largest trajectory induction residual is 1.00176e-11.
All seven final steps independently replay exactly from saved penultimate states.

The physical reconstructed total-energy changes, divided by the small analytic
wave energy, are -4.19458e-9, 3.85030e-8, -2.97506e-8 on the three meshes.
Corresponding grid changes are -2.71968e-9, 4.12707e-8, -2.99866e-8. They are
different measures and neither is claimed exactly invariant in finite precision.
The 60-step energy change is -1.10682e-6; additional steps eventually encounter
background/roundoff accumulation (480 steps: -2.56591e-7). The tenfold weaker
wave has relative physical energy change 3.29316e-6. Stable conservative-change
reconstruction makes quadrature 16/20 agree to approximately 1e-16 of wave
energy; the earlier less accurate background-subtraction audit is retained.

## Slow branches and precision limits

All branches are kept. At amplitude 0.1 m/s in equal-energy coordinates, the
slowest branch's B-rate error at 8/16/32 elements is 0.217142, 0.00374155,
0.00142150. The middle branch is 0.00115781, 0.000297148, 0.000296598.
The independent reference is linear; these are not all spatial errors.

At 32 elements and amplitudes 0.01/0.001, slowest-branch B errors are
1.52081e-4/5.61547e-5 and middle-branch errors 2.96613e-5/2.98077e-6.
The smallest-amplitude 16→32 element controls retain spatial convergence in
B/J. Density/pressure nonlinear residuals shrink with amplitude. However, the
weakest slow-branch velocity and E reach a numerical background-precision floor:
slowest-branch E at 32 elements rises from 1.38431e-6 to 1.38162e-5 as amplitude
falls tenfold. This limitation is retained, not called converged nonlinear
accuracy. Future weak-signal work should inspect constant potential gauge and
EOS difference evaluation before interpreting that floor physically.

## Stationary spectrum and strict refinement

The separately derived complete 360-DOF stationary matrix has spectral radius
2.47011516253e8 /s, maximum real/radius 1.82052e-15 and physical-energy symmetric
defect 1.51643e-13. Centered mixed nonlinear perturbations at 1e-4/1e-5 approach
it with energy-norm errors 8.29681e-7/8.29802e-9; full E errors are
1.89005e-6/1.89005e-8. Its RK4 imaginary-axis limit is 11.4506 ns for this
specific frozen mesh. It is not a general nonlinear or finer-mesh CFL rule.

The original full coordinate-difference attempt rejected a localized near-axis
toroidal magnetic coefficient (packed index 56, g=3.99965823e-6), with uniform
resting fluid. The initial pressure-only hypothesis was disproved by the saved
input; a frozen comment predates that diagnosis. Inner physical Faraday residual
4.23254e-11 passed, but the full coefficient backward residual 2.92848e-8 did
not. Both actual curl evaluations agree to 6.48e-19; no curl formula or physical
equation was changed. The original rejection replays exactly and is atomic.

New opt-in **RefinedRegularPlasma / RefinedKernelInduction** retain every
original full residual and its 1e-10 tolerance, tightening only the inner solve
when needed. The captured case passes at 9.56378e-13 using inner tolerance
1e-13. Complete coordinate matrices now finish at both amplitudes and replay
bitwise. Their independent linear-operator errors are 2.04814e-9/2.20696e-11;
physical-energy symmetry defects 2.89651e-9/3.12114e-11. Of 721 solves per
matrix, 16/62 require refinement, with up to 3/4 additional attempts. Every
original checked residual stays below 1e-10. The complete 120-step coarse wave
is bit-for-bit identical to the original arrays and non-time metrics. Exhaustion
still rejects atomically. The frozen original class remains available.

## Remaining gates

Continue with the opt-in refined assembly for a complete nonlinear contact
transit, positive-state/atomic-failure controls and actual time/mesh refinement.
That 100-microsecond transit is much longer than a fast-wave period; explicit
RK4 cost motivates a derived stiff time method, not skipped fast physics.
Close the original local physical enthalpy/Poynting interface budget; global
graph work and integrated energy are insufficient for that claim. Weak pressure
precision, earlier fluid accuracy limits, heated moving exact vacuum, physical
closures/ionization, actual coils/ends and converged formation remain open.
Neither these axisymmetric waves nor a future closed flux contour establishes
three-dimensional stability or particle confinement. Hardware remains abstract
prescribed waveforms, and the user's provisional lab envelope remains a set of
design variables, not a validated operating recommendation.
