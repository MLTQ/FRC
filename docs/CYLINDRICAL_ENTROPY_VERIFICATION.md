# Cylindrical joint transport: measured evidence

Final source14223ded4171ee70b3b790323913b1ff302d994df98a13c6328522273ebcf6ef
is frozen in outputs/cylindrical_entropy_verification. The full suite passes
391 tests in165.29s; the six focused tests also pass after the diagnostic
velocity-normalization correction. The independent audit covers17 original
2-D rate cases,6 radial spectra and17 unforced trajectories. No formation run
or production plasma adoption is included.

## Full original equations, including the axis

At fixed128 axial points, degree5 full acceleration relative L2 errors for
2/4/8/16/32 radial elements are5.93419e-4,1.81504e-5,5.68802e-7,1.79135e-8,
9.91917e-10. First-element radial errors are3.09383e-4,8.64907e-6,1.24677e-6,
9.05808e-8,5.87185e-9. The finest maximum **pointwise** relative radial error is
2.63e-7: the whole-domain norm must not hide this less favorable axis measure.
Degree3 also converges, more slowly near the axis. Radial DOFs are(p+1) per
element; these are not equivalent mesh counts to the previous face method.

Degree7 full error plateaus near8.18e-10 at128 axial points. Doubling axial
points reduces it to1.31660e-11/1.29162e-11 at16/32 radial elements, confirming
the axial truncation floor. Near-axis degree7 errors reach roughly1e-11 and
then become nonmonotonic; no unlimited-order or precision claim is made.

Grid physical H_t global defect is at most6.68e-15 relative to integrated
absolute local H_t. Density/angular/axial/entropy rate conservation is within
8.74e-15; independent local EOS/kinetic differentiation matches H-gradient
work within2.52e-16. Tests also check every radial element's boundary energy
balance. These identities alone are not accuracy evidence.

## Radial sound and unforced evolution

Analytic radial frequencies are c times the zeros of J1 divided by radius.
Degree3 fundamental relative errors at2/4/8 elements are4.04580e-5,6.92265e-7,
1.08595e-8; degree5 reaches9.10e-15. The largest dimensionless real eigenvalue
is3.02e-14 and physical acoustic energy-adjoint defect is below1.81e-15.
Full nonlinear RHS differences at amplitudes1e-3/1e-4/1e-5 confirm the
linearization, with a visible finite-difference roundoff limit at the smallest
amplitude. Matrix-exponential and RK4 time controls retain all raw matrices.

Three complete-transit contacts at16/32/64 axial points have continuous density
errors2.37128e-4,4.25734e-6,6.91715e-8. At fixed32 points, true grid-energy drift
is−1.50030e-11,−9.30478e-13,−5.79737e-14 for base/half/quarter step: factors
16.12/16.05. Reverse flow is retained separately, including its opposite signed
drift. Maximum conserved state-history drift across all17 runs is3.38e-15.

Rotating equilibria run for one radial sound-crossing time, not a full rotation.
At2/4/8 degree5 elements, density errors are1.06603e-8,3.10858e-11,5.42048e-12;
maximum spurious vr/c is4.25163e-9,8.39389e-10,3.83454e-11. Exact discrete
rotation balance is not asserted; the entropy-variable boundary is nonlinear.

Full-period Bessel acoustic runs at degree3 show density errors relative to
the initial perturbation2.59045e-3,3.62700e-4,7.45546e-5 at2/4/8 elements.
The reference is linearized, while the solver advances nonlinear equations.
At4 degree5 elements, reducing amplitude1e-4→1e-5→1e-6 reduces perturbation-
relative departure6.01946e-5→6.02175e-6→6.25897e-7. Absolute density error
therefore falls approximately quadratically, as expected for nonlinear
departure from a linear solution. Halving the time step at1e-4 changes that
relative error by only7.17e-12. Local entropy changes remain explicitly saved;
global entropy conservation is not a promise of exact local isentropy.

## Independent continuous integration and retained attempts

Every saved accepted snapshot is independently reconstructed with regular-axis
physical fields. Increasing quadrature32→48 radial and4x→8x axial changes
integrated energy by at most4.35e-16. The grid/continuous energy discrepancy is
larger: up to8.83e-11 over these histories. In the coarsest2-D manufactured
degree3 input it is2.53e-8; this too refines and is not hidden.

For rotating2/4/8-element cases, **continuous** energy drift is8.46457e-11,
2.21232e-13,1.11875e-14 despite grid-energy drift near roundoff. Continuous
angular-momentum drift is1.40247e-9,3.70942e-12,1.84170e-13. Maximum angular
grid/continuous difference falls1.46e-9→1.92e-11→2.14e-13. These are finite
representation errors, not physical sources or post-step corrections.

All records pass source/raw hashes, original-platform bitwise rate/spectrum/
last-step replay, and physical history/integration audits. The original attempt
under outputs/cylindrical_entropy_before_velocity_norm retains undefined
relative-velocity scores at zero-velocity turning points. The corrected score
uses RMS error/sound speed. All153 shared raw-array files replay bitwise after
that reporting-only change; neither the solver nor an evolved state was repaired.
All nine previous frozen source closures remain unchanged.

The six-panel figure is in outputs/cylindrical_entropy_diagnostics/report.png
(and PDF), with full saved-state continuous integrals in summary.json. Tests
and records are local. There is no configured Git remote or public publication.

## Remaining work

Verify fully mixed nonlinear r/z flow with temporal and mesh controls, preserve
the honest continuous/grid distinction, and derive the current/magnetic bridge
on this actual common state. Then verify full species/Maxwell/drift stress and
enthalpy/Poynting energy together. Heated free-edge exact vacuum, physical
collisions/conduction/ionization, real ends/coils, shocks/kinetic validity and
converged formation remain open. No axisymmetric equilibrium certifies3-D
stability or particle confinement.
