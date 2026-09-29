# Nonlinear compressible toroidal MHD verification

The bounded compressible core now evolves toroidal magnetic field, rotation and
conservative axial angular momentum, in addition to poloidal motion/field,
density and total energy. Scalar resistivity includes the new current components
and positive Ohmic work. **This is a verified component milestone; nonlinear Hall
coupling and the physical plasma/vacuum formation interface remain incomplete.**
The original low-pressure argon target is unchanged.

The full suite passes **224 tests in 85.76 s**, including 16 new checks. They
cover exact annular kinetic energy/angular momentum, rotation/pinch balance,
independent unbalanced force, both-order reduction to the old poloidal RHS,
nonlinear transport, torsional propagation, conducting-Bessel decay, compatible
positive resistive work, regular-axis polynomial diffusion, failure atomicity,
nondefault gamma and independent physical-volume reconstruction error.

## Complete reproducible study

All **33 cases** completed: nine rotation/pinch equilibria, eight contact/wave
mesh refinements, three nonlinear contact controls, four wave controls, five
radial diffusion meshes through nr=128, and four time caps. Full configurations,
lossless initial/final/reference states, initial RHS and accepted-step histories
are in `outputs/toroidal_verification`. Numeric-column histories retain every
sample while avoiding repeated object overhead. Sources/configurations/archive
identities and complete terminal histories are verified. Frozen study SHA256:

`24ce9b8cdb121167d78ab2ee75792bb9ace362695dd852cbc2384065f4dc69f2`

Report/render sources and figures have separately recorded hashes. Historical
Hall, current and moving-boundary studies remain intact. Reproduce from the
repository root:

```sh
MPLCONFIGDIR=/private/tmp/plasma-bench-mpl XDG_CACHE_HOME=/private/tmp/plasma-bench-cache .venv/bin/python -m examples.toroidal_measurements
MPLCONFIGDIR=/private/tmp/plasma-bench-mpl XDG_CACHE_HOME=/private/tmp/plasma-bench-cache .venv/bin/python -m examples.verify_toroidal
.venv/bin/python -m pytest -q
```

All fixtures are explicitly SI mathematical references with prescribed ideal-gas
and scalar-eta parameters. They do not select laboratory operating parameters
or bypass the original regime/edge uncertainties. Equations and numerical
stabilization are specified in [TOROIDAL_MHD_NUMERICS](TOROIDAL_MHD_NUMERICS.md).

## Nonlinear transport, stresses and waves

Independent rigid-rotation, linear-pinch and combined equilibria retain their
pressure/centrifugal/magnetic balance. Across nr=8/16/32, maximum spurious radial
volume-RMS velocity after 0.2 s is **3.99221e-16 m/s**. An unbalanced reference
checks the exact acceleration `(Omega²-2c²)r`, including axis and wall. This
special polynomial preservation is not a claim of arbitrary equilibrium balance.

The finite-amplitude rotating entropy wave has rho=1+0.5 cos(pi(z-Vt)),
Omega=0.7/sqrt(rho), V=0.8 m/s, Btheta=sqrt(mu0)*0.4r, Br=Bz=0 and fixed
quadratic p(r). True continuum conserved averages initialize it; no energy
adjustment forces numerical pressure to agree. The solution translates through
the periodic domain and exercises nonlinear angular-momentum/energy transport.

Two different norms are retained: errors of cell averages/recovered regular
coefficients, and independent physical-volume errors of the actual limited
reconstruction against exact continuum point fields. The latter includes
subcell error. Paired tensor Gauss orders 6/10 change any reported physical
relative norm by at most **8.21825e-17** in these smooth cases.

| Axial cells | Contact reconstructed density error | Contact reconstructed rotation error | Torsional reconstructed rotation error |
|---|---:|---:|---:|
| 24 | 0.00620541 | 0.00848995 | 0.0197280 |
| 48 | 0.00165455 | 0.00188532 | 0.00543534 |
| 96 | 0.000449627 | 0.000446635 | 0.00158995 |
| 192 | 0.000126545 | 0.000125193 | 0.000462097 |

These are relative volume L2 norms at t=0.4 s. Contact density/rotation cell-
average errors at nz=192 are 0.000123532/0.000122588, explicitly distinct from
the reconstructed-field errors. Its numerical radial velocity is 3.95150e-6 m/s
RMS, decreasing with refinement. Doubling radial cells at nz=96 gives little
change, as expected for this radially simple exact state.

A stronger 0.8 density contrast remains positive but has larger nz=96 rotation
cell error 0.00306026; it is not declared converged from one mesh. First-order
nz=96 contact density/rotation errors are 0.0300126/0.0319303, exposing the
explicit numerical diffusion. No density/pressure floor or hidden fallback is
used. Minimum accepted final pressure over the protocol is 0.757690 Pa.

The torsional wave is an independent linear Alfven reference, with regular
vtheta and Btheta proportional to r. Its relative nonlinear correction is
small at amplitude 1e-4 rad/s. Amplitudes 1e-3/1e-4/1e-5 give reconstructed
density departures 1.13299e-7/1.13299e-9/1.13300e-11, as expected for quadratic
feedback. The rotation discretization error is essentially unchanged, separating
numerical error from that physical nonlinear response. A third axial harmonic
is also recorded; first-order nz=96 wave error is 0.0516970 in cell values.
Neither a linear wave nor this smooth nonlinear contact verifies strong shocks.

## Toroidal diffusion and energy

The independent conducting-wall mode is Btheta=C J1(qr) exp(-Dq²t), where
J0(qR)=0. C=1e-7 T, D=0.1 m²/s and t=0.2 s make fluid feedback asymptotically
small. The radial current mass is a compatible positive tridiagonal operator.
It recovers the exact 8Db diffusion of g=a+b r² including the axis; diagonal
lumping would give the wrong axis coefficient. No such lumping or energy
correction is used. Semidiscrete magnetic work equals positive resistive heat;
SSPRK2 time error and actual continuum accuracy remain separate checks.

| Radial cells | Toroidal field cell error | Reconstructed field error | Axis g error magnitude | Physical heat error |
|---|---:|---:|---:|---:|
| 8 | 7.09646e-4 | 1.05666e-3 | 4.56566e-4 | 5.34243e-3 |
| 16 | 1.45390e-4 | 2.29640e-4 | 9.19367e-5 | 1.33219e-3 |
| 32 | 3.30830e-5 | 5.49162e-5 | 2.03664e-5 | 3.32373e-4 |
| 64 | 7.94882e-6 | 1.35231e-5 | 4.77566e-6 | 8.29895e-5 |
| 128 | 1.95079e-6 | 3.35910e-6 | 1.15482e-6 | 2.07328e-5 |

All are relative errors; the axis quantity is the first annulus's averaged g,
not a falsely nonzero Btheta at r=0. Physical heat is compared with the exact
continuum magnetic-energy loss. Toroidal flux may change through finite Ez on
the coordinate axis under resistivity; the axis is not an artificial conductor.

Across the full study, final mass and total-energy relative changes are at most
4.52346e-15 and 1.13086e-15. Relative angular momentum change is at most
2.05361e-16 for cases with a nonzero resolvable total (not normalized by tiny
roundoff totals of zero-mean waves). These budgets do not establish field accuracy.
The independent time-cap comparison uses the same nr=8,nz=32 nonlinear contact
at t=0.05 s and a dt=6.25e-5 reference. Density halving ratios are 4.1723/3.9905;
angular-momentum ratios are 3.7153/4.3801. Limiter switching means they are not
claimed to be an exact asymptotic factor four in every variable.

![Toroidal accuracy](../outputs/toroidal_verification/toroidal_accuracy.png)

![Actual nonlinear transport](../outputs/toroidal_verification/rotating_transport.png)

## Remaining integration

Next implement [NONLINEAR_HALL_PLAN](NONLINEAR_HALL_PLAN.md): compatible Hall
induction/Poynting work on the new compressible state, independently checked
finite-sound-speed cylindrical modes and nonlinear transfer. The older
incompressible wave reference is not an exact compressible oracle. Electron
pressure and finite-mass current/thermal transport require the full consistent
parent-species reduction; they remain incomplete.

The thermal free-edge conflict, spatial coil work, actual vacuum/chamber ends,
shock/remapping robustness, low-density kinetic ordering, converged formation
and experimental UI remain active requirements. No new data demonstrates FRC
formation, 3D stability or particle confinement. Stage 2 remains guarded.
