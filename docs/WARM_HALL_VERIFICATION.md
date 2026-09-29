# Warm Hall integration: measured progress with explicit failures

The full production regression suite passes **293 tests (94.10 s)**.
All 58 main cases finish and pass source, configuration, archive, history and
independent physical-energy/field audits. Five supplemental cases finish; the
sixth fails at finer radial resolution. Three dedicated replays preserve that
failure and its time-step controls. Thus there are **63 completed trajectories/
spatial measurements and three rejected trajectory records**, not 66 successes.
The new pressure-only principal-symbol counterexample restricts this split
model's use. This is verified research progress, not FRC formation readiness.

## Scope and reproducibility

The model evolves moving ion density/momentum, poloidal/toroidal fields and
separate equal-gamma species temperatures, using conservative Hall/pressure
Poynting power and electron-current enthalpy flux. The cylinder is periodic in
z, with fixed psi and Btheta(R)=0, free Ez(R), and no normal fluid/charge flux.
It is not a conducting electric wall or a true-vacuum boundary. Coupled eta=0;
resistivity appears only in the separate boundary diffusion test.
WARM_HALL_NUMERICS states all reductions and local numerical heat limitations.

Frozen engine/reference/measurement fingerprint:

`ff529e0fab60fb0888806589056142c7040b3c66fcff0cfde669635e1fbb7c20`

The 51 transitive source files, Python/NumPy/SciPy versions, full cases, accepted
histories and lossless field SHA256 values are recorded. Reports separately
hash their audit/plotting code and three visually inspected figures. Results
are in `outputs/warm_hall_verification`, `outputs/warm_hall_controls`, and
`outputs/warm_pressure_failure`. Failure records retain false trajectory
completion flags and separately labeled rejected arrays.

The previous 36 moving-ideal records retain unchanged engine identity
0e37c39a… . The prior 46 electron-pressure records now use historical source
verified against commit 419c951, since the wall-density guard changed. All 63
older cold Hall records were source-verified against 16969a9 after changes to
shared mass/Poynting interfaces. These records were not relabeled as new runs.

## Local pressure work and boundary controls

The independent nonlinear reference supplies analytic pressure work, induction,
thermal rate and physical face flux integrals. Their continuum energy identity
agrees within 1.15e-13 relative to the pressure-work scale. Gauss 12→18 changes
the reference cell energy rate below 7e-16 in that normalization.

| nr | Conservative energy-rate relative L2 | Local ion-energy defect / pressure work | First-annulus energy error / pressure work |
|---:|---:|---:|---:|
| 8 | .0940091 | .108611 | .0836177 |
| 16 | .0248245 | .0284782 | .0215521 |
| 32 | .00629279 | .00720739 | .00542908 |
| 64 | .00157868 | .00180742 | .00135985 |
| 128 | .000395013 | .000452206 | .000340122 |

The first column divides by the exact energy-rate norm; the other two divide
by pressure-work norms. Gamma1.4 and uniform-density controls also refine.
For the barotropic reference the exact energy rate is zero, so relative-to-zero
errors are null. Its residual divided by the nonzero pressure-work norm is
.0169391→.00428555→.00107459 at nr16/32/64. Curl-free induction alone does
not guarantee exact discrete enthalpy/Poynting cancellation.

The new full-Ez boundary satisfies C grad=0, its weak work matches the dual
current at every DOF, and constant-scalar current remapping preserves charge.
Varying wall n/pe permits nonzero Ez with zero wall magnetic/enthalpy power.

The separate zero-J1 toroidal resistive mode has a **nonconvergent instantaneous
first-annulus derivative error**: signed error approaches -.152131 at nr128.
Volume g-rate errors decrease because the faulty region shrinks. After decay
exponent .02, evolved g errors at nr16/32/64/128 decrease toward 5.96195e-6,
but the two-half-step control changes the finest error to 3.66451e-6. Magnetic-
weighted coefficient errors are 1.27723e-5 and 7.64079e-6 respectively. These
are coefficient norms, not pointwise B accuracy or an erased axis defect.
Positive resistive heat agrees with physical magnetic energy loss. Coupled
warm resistivity stays disabled.

## Actual warm waves

Uniform background: R=.035 m, L=.24 m, n0=1e20 m^-3, argon ion mass,
B0=.005 T, Te=5 eV, Ti=1 eV, gamma=5/3, reference velocity amplitude 1 m/s.
Independent Bessel eigenmodes satisfy continuity, momentum, Faraday and both
species thermal equations. Each evolves through phase .5 rad. All errors below
are physical L2 norms normalized by the perturbation, not the large background.

| Branch | nr | Electron pressure | Ion pressure | Magnetic field |
|---:|---:|---:|---:|---:|
| 0 (slow) | 4 | .801711 | 2.47144 | .970361 |
| 0 | 8 | 1.10769 | 6.73716 | .493712 |
| 0 | 16 | .902996 | 4.44172 | .116057 |
| 0 | 32 | .164298 | .829787 | .0407667 |
| 1 | 4 | .197019 | .927567 | .328317 |
| 1 | 8 | .0537335 | .262316 | .165546 |
| 1 | 16 | .0136194 | .0654738 | .0829770 |
| 1 | 32 | .00339890 | .0165829 | .0415157 |
| 2 | 4 | .193546 | .999840 | .312328 |
| 2 | 8 | .0537838 | .278136 | .157426 |
| 2 | 16 | .0138332 | .0692166 | .0789185 |
| 2 | 32 | .00347082 | .0173594 | .0394869 |

Fast species pressure errors approach second order; physical magnetic
representation remains approximately first order. Finest branch2 electron/ion
temperature errors are .00803288/.0430700. Slow branch0 remains insufficiently
accurate: even nr32 has 16.43% electron-pressure and 82.98% ion-pressure error,
with ion-temperature error 234.86% of its small perturbation. Its electron
pressure RMS error is 1.89717e-5 Pa. Small absolute values cannot justify hiding
the relative mode error. Finest slow evolution takes 4,344 accepted fluid steps.

At nr16, half CFL changes slow pe RMS by 1.36107e-6 Pa; its exact-reference
relative error remains .900427 versus .902996. Reducing amplitude tenfold
retains relative error .903035; normalized nodal pressure differences are
.00204560 of the baseline perturbation. These controls do not explain away
the large spatial error. Opposite axial signs, near-cold electrons and fast
amplitudes .1/10/100 m/s are retained separately. Independent physical norms
at orders 8/12 differ by at most 1.62e-11; electric-quadrature versus independent
norms differ below 2.58e-10 in relative-L2 units.

## Time, nonlinear evolution and conservation

Fixed-grid pe difference ratios under step halving are 4.00782/4.04662 for
isolated pressure, 3.99940/4.05812 for the fast wave, and 4.00586/4.23940 for
the nonlinear coupled fixture. Magnetic ratios are separately retained and
range about 3.82–4.09. This verifies smooth temporal convergence in these
windows, not general stability of the pressure split.

The nonlinear variable-density, variable-wall-pressure fixture reaches 1 µs
at nr8/16/32. Continuous pe differences from nr32 decrease .373756→.0836107 Pa,
or .187270→.0418930 of the fine trajectory's pressure change. Physical magnetic
differences are .222484→.0998723 of the fine field change. There is no exact
nonlinear time oracle; these are nested-mesh comparisons including magnetic
representation error, not an asserted converged formation solution.

Maximum main accepted-history relative drifts are 5.41e-16 mass,
2.15e-15 total energy and 5.82e-15 linear electron entropy proxy. Physical
logarithmic entropy changes are tracked separately: up to 1.33542e-6 per
particle in the coarse nonlinear time control. The earlier bulk entropy defect
and lack of a shock entropy inequality remain. Global conservation does not
bound local ion-energy errors or establish species accuracy.

## Failed fine pressure control and next step

In the stationary-ion pressure subproblem, RMS spurious ion-pressure changes
at nr8/16/32/64 are .0607924/.0163608/.00416953/.00116477 Pa. However, maximum
errors are .461809/.130961/.0350927/.0920654 Pa: the finest maximum worsens.
Finer time steps retain this. At nr128, ion pressure becomes negative at about
.95 µs in three time-resolution controls, despite positive electron pressure
and conserved energy. These trials are rejected and archived, not repaired.

WARM_PRESSURE_SPLIT_LIMITATION derives and preserves the independent negative
principal discriminant of the isolated pressure operator and all rejected
states. The next gate is a joint fixed-density Hall/pressure integrator with
consistent magnetic midpoint, thermal mean, local energy flux and atomic
validation. The observed axis failure also requires independent axis-local
analysis; the continuum symbol is not proof of its particular cause.

Slow-wave accuracy, resistive-axis consistency, physical collision/heating,
finite-mass current transport, conduction/ionization, original heated edge,
true vacuum, real gas ends, spatial coils, shock/remapping robustness, kinetic
ordering and converged FRC formation remain unfinished. Stage2 stays guarded
and the full goal remains active.

Reproduce the main measurements and report with:

```sh
.venv/bin/python -m examples.warm_hall_measurements
.venv/bin/python -m examples.verify_warm_hall
```

For the follow-up evidence, run both `warm_hall_controls --group slow` and
`warm_hall_controls --group pressure` (the latter deliberately reproduces the
reported exception), then `warm_pressure_failure`, `pressure_split_symbol`
and `verify_warm_hall_controls`, each as a Python module. Do not interpret the
control driver's nonzero exit as a completed trajectory or skip its replay.
