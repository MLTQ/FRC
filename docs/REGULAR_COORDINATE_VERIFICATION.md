# Regular-coordinate fluid: measured verification

This is a smooth, positively filled, fixed-radius, periodic fluid component.
It has no magnetic/current coupling, moving vacuum boundary, shock treatment or
physical formation result. Equations, axis degrees of freedom and boundaries
are in REGULAR_COORDINATE_NUMERICS.md. All results below use saved states and
their declared physical reconstruction, without filtering or energy repair.

## Why the representation changed

The previous Gauss-in-r source conserved its grid energy but used a different
regular polynomial in r² to define continuous physical fields. Mixed unforced
radial/axial/rotating flow exposed large between-point velocities. The old
14-case study is preserved in outputs/mixed_euler_verification under source
a3d0b0c387bd4696fca958c20038d41ad7725746f4809d4052f80c36493e07b6.
All trajectories completed and replayed; this does not mean their coarse
physical representation was accurate.

At degree 5 and 2/4/8/16 radial elements (16/32/64/128 axial points), final
old grid/continuous energy gaps were 0.0667539, 0.0253020, 0.000119298 and
8.41354e-7. Maximum sampled velocity component divided by sound speed was
8.9899, 14.3621, 2.39269 and 0.296464, despite point maxima near 0.295.
Time refinement did not remove the representation error.

An independent polynomial integration in outputs/axis_interpolation_diagnostics
compares the actual first-element kinetic Gram matrix with the grid metric.
At degree 5, the largest generalized eigenvalue is 576461.02; the ratio is
invariant under element scaling. Applying the physical work adjoint of f(r)=r
and then the old lumped inverse gives relative force error 0.66543. This rules
out that naive coupling bridge, not every possible Gauss-r formulation.

The new active first-element velocity polynomial has degree p−1 in x=r²/2.
Its constant-density physical kinetic metric is exactly the positive nodal
metric: degree 3/5/7 eigenvalues differ from one by at most about 3e-14 and
the same force errors are at most 5.26e-15. All three independent metric/force
records have verified hashes and bitwise replayed arrays and metrics.

## Frozen new evidence and audit scope

New source identity:
06dfff4762b721a154d9b8522f935427da0a13f335ac74a69549a08cdc958e7e.
The directory outputs/regular_euler_verification contains 8 original-equation
rate cases, 8 full acoustic spectra and 14 mixed evolution cases: 30 audited
records. Together with the 14 old trajectories, all 44 records are retained.
The full regression suite passes **398 tests in 181.18 s**.

Audits verify source/raw hashes, original-platform bitwise rates and spectra,
initial/final evolution rates and the last accepted state/work increment.
Every saved continuous energy/angular history is independently reconstructed
from active degrees of freedom. The accumulated manufactured source-work
ledger before the final step is retained but not independently replayed in
full. Audit success establishes replay and measured checks, not blanket
accuracy of every mesh or acoustic branch.

Two physical quadratures (32 radial sampling elements with 24/40 Gauss points
and 256/512 axial samples) agree in new energy within 2.76914e-16 across saved
snapshots. Independent local physical energy derivatives agree within
3.32714e-16; conserved-source balance remainder is at most 6.61642e-15.
Ten earlier source closures, the old mixed comparison and the new source are
verified unchanged in regular_euler_diagnostics/source_preservation.json.

## Original equations and acoustic modes

Independent rate fields use the original regular-coordinate Euler equations,
not the solver flux identity. At degree 3, first-element active q-rate errors
are 1.64243e-5, 2.44402e-7, 3.77987e-9 and 7.25635e-11 at 2/4/8/16 elements.
Degree 5 reaches 6.32753e-10 and then the fixed axial-grid floor near 4.21e-11.
Maximum global grid physical energy-rate defect is 2.02609e-15.

The independent linear operator uses delta_t=−2 D_x(x a),
a_t=−R² D_x delta with a=R*u/c and a positive active acoustic metric.
At degree 5 and four elements, its first three Bessel-frequency relative
errors are −4.869e-13, −1.202e-8 and −5.971e-7; eight elements reduce them to
approximately −3e-15, −1.209e-11 and −6.133e-10. Actual nonlinear RHS
differences perturb both real and imaginary mode components. Their error
decreases from 3.34566e-8 to 3.34691e-10 before roundoff at 1e-5 amplitude.
RK4 errors against the matrix exponential over two periods decrease by about
16 under step halving.

**The degree-3 third sorted positive mode is not the third Bessel root.**
Its roughly 29–30% mismatch persists under refinement. Full spectra retain
this extra branch; only verified mode claims above are supported. No blanket
spectral accuracy or production adoption is claimed.

## Matched mixed nonlinear evolution

Both methods start from the same analytic density, all three velocities and
two species pressures, for one radial sound-crossing time R/c. The independent
time-dependent exact field supplies manufactured forcing; its wall radial
velocity is zero and p_r=M*n*vtheta²/r. Unforced cases use the same initial
field but have no known exact final solution. Base time steps differ between
the schemes; each has independent half/quarter-step controls.

| Radial elements / axial points | New manufactured density relative L2 | Central velocity RMS/c | Unforced final grid/continuous energy gap |
| --- | ---: | ---: | ---: |
| 2 / 16 | 1.82947e-5 | 9.69768e-6 | 5.72337e-8 |
| 4 / 32 | 3.62103e-7 | 1.69433e-7 | 5.16148e-9 |
| 8 / 64 | 7.05634e-9 | 2.72857e-9 | 3.63098e-10 |
| 16 / 128 | 1.64671e-10 | 4.29509e-11 | 1.70248e-11 |

The central norm covers the fixed physical region r<R/8. Unforced maximum
sampled component speed/c stays near 0.29646 on every mesh. New grid/continuous
energy differences throughout all histories are at most 9.04420e-8; these are
measured discrepancies, not exact continuous conservation.

Unforced continuous energy drift refines from −5.72815e-8 to −1.70396e-11;
continuous angular drift from 1.06029e-8 to 4.85773e-14. At four elements,
grid energy drift is −3.06793e-12, −1.38457e-13 and −7.19977e-15 for base,
half and quarter steps. Energy-drift ratios are not exactly 16. Actual field
time differences have ratios 15.88 (density) and 16.09 (velocity), consistent
with fourth-order time convergence in this control.

**Free-flow spatial convergence is only roughly second order in this study.**
Successive common-grid velocity RMS/c differences are 7.06396e-4, 1.85237e-4
and 4.70508e-5; density differences are 5.82540e-4, 1.63202e-4 and 4.09450e-5.
The finer numerical result is not an exact solution. Higher wall compatibility
or evolved regularity may matter, but that explanation is not established.
The stronger manufactured convergence does not erase this limitation.

The inspected six-panel report and complete field/representation/time
comparisons are in outputs/regular_euler_diagnostics. No data were filtered.
The next gate is an independently verified variable-density current/magnetic
work and torque bridge on this actual state. Full species/drift/Maxwell stress,
enthalpy/Poynting transport, physical closure, heated moving exact-vacuum edge
and converged FRC formation remain open.
