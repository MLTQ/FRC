# Electron pressure and entropy verification

**266 tests pass (88.54 s).** All 46 declared records are complete: 17 spatial
cases, 28 successful evolutions and one deliberately strict nonlinear-solver
failure. The new component preserves conservative charge flux, couples electron
thermal drift to pressure induction, and accounts for discrete physical entropy
alongside magnetic/electron energy. It does not yet evolve warm ions or the full
compressible Hall system. ELECTRON_PRESSURE_NUMERICS states the equations,
changed scalar placement, trace boundary and remaining physical scope.

## Independent spatial and boundary checks

Analytic positive density and temperature have nonparallel gradients; an
independent smooth toroidal field supplies current and entropy advection. Both
pressure-driven induction and thermal transport are differentiated from physical
continuum fields. Exact current face integrals use differences of r Btheta;
they are not computed with the numerical current map. Physical annular averages
at orders10/14 change the induction reference by at most 6.79e-16 relative.

| nr (nz=2nr) | Battery gdot L2 error | First-annulus error | Electron thermal-rate error | Axial current face-integral error |
|---:|---:|---:|---:|---:|
| 8 | 0.0252946 | 0.0249413 | 0.0985281 | 0.0258781 |
| 16 | 0.00636035 | 0.00626256 | 0.0258587 | 0.00711227 |
| 32 | 0.00159239 | 0.00156732 | 0.00654621 | 0.00188419 |
| 64 | 0.000398242 | 0.000391935 | 0.00164173 | 0.000485857 |
| 128 | 0.0000995696 | 0.0000979904 | 0.000410757 | 0.000123411 |

Radial current face-integral error reaches 6.50492e-5, and first-axis axial
current error reaches 2.50986e-5. Stronger density contrast and gamma=1.4 retain
the observed second-order rates. These integrated-current and coefficient norms
are not interchangeable with continuum pointwise magnetic-field accuracy.

Uniform-density pressure has exactly zero induction. Isothermal variable-density
pressure has only roundoff-level residual: absolute gdot RMS errors are
2.51911e-10, 6.58314e-10 and 1.31973e-9 T/(m s) at nr16/32/64. No relative
error is invented for a zero reference. A separate adiabatic barotrope also
has zero curl to roundoff. The growth of this tiny absolute residual under
differentiation is retained; it is not physical field generation.

Independent polynomial current tests cover the wall trace. Btheta=C r yields
uniform Jz=2C/mu0 all the way to the wall. Omitting its nonzero wall trace adds
a spurious return-current sheet and gives the wrong sign in the boundary dual
cell. The general dual work test retains nonzero prescribed wall current and
Poynting power, rather than falsely declaring every wall thermally closed.

## Actual electron heat translation

Uniform n=1e20 m^-3 and Btheta=0.01 r T/m give constant axial electron velocity
ve=-2C/(mu0 n e), while bulk velocity is absent. A 30% pressure wave with
mean 80 Pa travels through one full 0.24 m axial period. Its pressure electric
field is curl-free; Btheta remains exactly unchanged in the saved runs.

| Axial cells | Nodal pressure perturbation error | Continuous reconstructed pressure error |
|---:|---:|---:|
| 16 | 0.168907 | 0.167636 |
| 32 | 0.0425101 | 0.0425302 |
| 64 | 0.0106488 | 0.0106771 |
| 128 | 0.00266371 | 0.00267227 |

The independent physical norm integrates the continuous Q1 pressure field
between nodes; order 8/12 changes this norm by at most 6.16e-17. The finest error
is **0.2672% of the pressure perturbation**, not of the large pressure background.
Opposite current/propagation sign agrees to rounding precision. The 70%
amplitude stress at nz64 has physical pressure error 0.0123249; positive
pressure alone is not used as an accuracy certificate.

Fixed-grid time controls at nz32 compare 64/128/256 steps to 1024. Pressure
differences divided by the fine-reference **change** are
0.122013/0.0301758/0.00718658, giving halving ratios **4.04340/4.19891**. Because
the analytic full-period solution returns to its initial state, this particular
denominator is mostly the reference's remaining spatial error; the physical
analytic errors above are reported separately. The magnetic time denominator
is zero and correctly recorded as null.

![Spatial, physical transport and time errors](../outputs/electron_pressure_verification/convergence.png)

## Nonlinear pressure battery and thermal feedback

The battery begins at Btheta=0 with smooth electron temperature near 5 eV and
nonuniform fixed positive density. Misaligned density/temperature gradients
create toroidal field; the resulting poloidal current advects electron entropy
and changes its pressure. No external magnetic forcing or thermal repair is
inserted. This is evolution of the **isolated pressure/entropy operator**, not
full plasma formation or a solution with the omitted Hall and bulk forces.

After 1 microsecond, magnetic energy at nr8/16/32/64 is
19.5070/19.8952/19.9889/20.0123 microjoules. Electron thermal energy decreases
by the same amount, within the measured global budget. At nr64 minimum electron
pressure is 49.8269 Pa. Nonlinear initial finite-Btheta and stronger-density/
temperature controls also complete. The strong nr64 case generates 278.219
microjoules with minimum electron pressure 39.9332 Pa.

Against the restricted nr64 reference, nr8/16/32 battery magnetic differences
are 0.0233811/0.00558931/0.00112172 of the reference change; dual pressure
differences are 0.101808/0.0258423/0.00528923. With nonzero initial magnetic
field, the corresponding nr32 differences are 0.00111263/0.00507041.
g is restricted by exact annular volumes; pressure is sampled at coincident
dual nodes. These are declared discrete-state comparisons, not proof that nr64
is an exact continuum solution or that its between-cell Btheta is higher order.

Fixed-grid nonlinear time steps 10/20/40 against320 yield pressure errors
4.58350e-4/1.14369e-4/2.82631e-5 of the reference change: ratios
**4.00765/4.04658**. An independent adaptive DOP853 test verifies temporal
order against the same semidiscrete ODE. Changing the arbitrary proxy density
scale from 1e20 to 7e18 changes evolved magnetic/pressure results by only
1.0030e-14/5.8183e-14 of their reference changes.

![Actual pressure battery and energy exchange](../outputs/electron_pressure_verification/battery.png)

## Conservation, finite precision and audit

Maximum spatial backward errors are **1.38766e-16** for local charge,
**1.80985e-16** for magnetic/thermal work, and **1.74751e-16** for physical
logarithmic entropy rate. Denominators use actual absolute summation scales,
not near-zero net exchange. Raw extensive entropy-rate roundoff can look large
because particle counts are about 1e17; both normalization and units are explicit.

Across all accepted histories, maximum energy drift is **5.70516e-15** relative;
the conserved linear entropy proxy changes at most3.75632e-16 relative. Maximum
discrete physical logarithmic entropy drift is **3.35616e-12 per particle**.
The last includes accumulated nonlinear solve tolerance over long drift runs;
it is not described as exact floating-point conservation of a reconstructed
continuum integral. Initial/final electron energy and entropy are independently
replayed from physical pressure and density, not just copied solver diagnostics.

The default nonlinear tolerance is 1e-10. At nr16, changing it to 1e-8 changes
magnetic/pressure trajectories by 1.05226e-10/5.46422e-10 of their reference
changes;1e-12 changes them by 2.79424e-13/2.22662e-12. The deliberately strict
nr64 first battery step at tolerance 2e-13 fails with residual9.95419e-13 after 80
iterations, retaining its exact initial state and zero accepted time. Larger
iteration-budget checks are separately recorded. The observed last-bit thermal
oscillation is amplified by differentiated axis coefficients; demanding a
residual below this observed precision level does not improve the solution.
The accepted lower-resolution/finer-time comparisons remain the accuracy evidence.

Source/configuration/individual/collection/archive identity, all finite arrays,
true accepted times, pressure positivity, solver residuals and volume conventions
pass audit. Both figures were visually inspected. Frozen engine/reference/
measurement SHA256:
`5cccf843cfe6ebbf6005d6b87261db7e9ae9decc1277a624b2ec90de8bbc78a7`.
Report source hashes, supplementary iteration controls and figure hashes are
separate in `outputs/electron_pressure_verification/report.json`. Reproduce with:

```sh
.venv/bin/python -m examples.electron_measurements
.venv/bin/python -m examples.verify_electron_pressure
```

All jobs are terminal. The next step is conservative primal/dual mass and
electron-energy coupling to moving ions, with compatible total pressure and
boundary traces. Collision exchange, resistive heating partition, thermal
conduction/ionization, finite electron inertia, physical vacuum/edge/ends,
spatial coils, shock/remapping robustness, converged FRC formation and the
dynamic experiment interface remain required. This is verified progress,
not completion of the active objective.
