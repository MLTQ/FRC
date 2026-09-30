# Joint transport verification checkpoint

The full suite passes **385 tests in 164.90 s**. The new point-state transport
component is separate from the running common-radial-face prototype and the
coupled plasma solver. Equations, the local physical-energy flux proof, state
meaning and uncompleted cylindrical/current integration are in
TWO_ENTROPY_NUMERICS. Source identity:
`924a586a84df69c28fe1493bf73e285efcecc02e63b6825926f5be232f94170d`.

All **12 original-equation rate cases and 13 unforced trajectories** pass source/
raw-array and bitwise replay audits. Independent physical sampling integrates
the declared continuous trigonometric state at 8x and 16x resolution. The maximum
energy difference between those two integration resolutions is 1.92e-16.
Maximum nodal-versus-continuous energy discrepancy among saved states is
6.70e-14; this separate error is retained rather than presumed zero.

## Original equations and conservation

Three arbitrary velocity components, nonuniform n and independently varying
electron/ion pressures are included in the instantaneous manufactured reference.
Raw rate keys radial_momentum/azimuthal_momentum label its two transverse
Cartesian components; this axial reference includes no cylindrical curvature.
Orders2/4/6 and 16/32/64/128 points are tested. The order6 physical-energy-rate
relative errors are 9.97242e-4, 2.16930e-5, 3.74766e-7 and 6.01370e-9.
Finest separate electron/ion pressure-rate errors are 4.01490e-10/1.19765e-8;
temperature-rate errors are 1.25633e-9/2.85462e-8. The ion errors are larger
and are not hidden in a combined species norm.

Across all rate cases, normalized global six-state conservation error is at
most 1.98e-16. Independently differentiated physical H_t sums to at most
1.02e-16 relative to its absolute local rate. Its agreement with the local
conservative energy-flux divergence is within 2.13e-14. These are spatial
identities in the declared integration norm, not claims of time accuracy.
Uniform-velocity contact rates have no added force from an energy correction.

## Full unforced trajectories

The original full-transit two-temperature contact retains unequal gammas and
the original physical n, pressures and M=mi+me. The nonlinear simple wave uses
its original equal-gamma analytic reference and duration. No analytic forcing,
energy reset, damping or density/temperature floor is applied.

| Axial points | Contact continuous density error | Wave continuous density error |
|---:|---:|---:|
| 16 | 2.37122e-4 | 5.14477e-5 |
| 32 | 4.25714e-6 | 1.12034e-6 |
| 64 | 6.92316e-8 | 1.92745e-8 |
| 128 | 1.09227e-9 | 3.15956e-10 |

These norms divide by the full physical density; separate density-perturbation
norms are retained in every record (finest contact/wave: 5.26361e-9/4.48378e-9).
Axial integral histories are per unit cross-sectional area. Finest contact nodal electron/ion pressure
errors are 2.10068e-9/3.10001e-9 and temperature errors 2.20153e-9/3.76193e-9.
Finest wave pressure errors are 5.22047e-10 for both species; corresponding
temperature errors are 2.08943e-10. Continuous combined-species norms, initial
representation and separate species endpoint norms are all available.

Maximum normalized drift of density, all three Cartesian momenta and both
entropy densities over all accepted histories is 1.79e-15. This is not yet
a cylindrical angular-momentum claim. All 13 trajectories complete with positive
physical nodal density and species pressures; independent sampled density is
also positive. No shock or global positivity theorem follows from these tests.

Global entropy conservation is distinct from exact local isentropy. Although
the analytic simple wave has constant specific entropy, the flux produces
small local variations. Maximum electron/ion specific-log-entropy errors are
about 1.15140e-6, 3.07600e-8, 6.10517e-10 and 9.96e-12 at 16/32/64/128 points.
Their rapid refinement is verified, but exact pointwise preservation is not
claimed. This is a retained accuracy limitation for future cylindrical work.

## Remaining energy error is temporal in this axial component

At 32 points, the full-transit contact's nodal physical-energy drift is
−2.46222e-11. Half and quarter steps give −1.52514e-12 and −9.40984e-14:
improvement factors 16.14 and 16.21. Independent continuous-energy drift tracks
these values within the retained quadrature discrepancy. Finest 128-point
contact/wave nodal energy drifts are −9.52407e-14/−8.08574e-14.

Reverse-contact energy drift is positive, 2.37207e-11, while the reversed wave
matches its forward counterpart. These reversible truncation effects are
retained; they are not physical heating. The contact density errors barely
change under fixed-grid time refinement, so its remaining field error is
spatial even though the energy drift is mainly temporal.

This addresses the earlier axial translation obstruction by changing the joint
mass/momentum/entropy transport. It does not retroactively give the old mixed
cell/face method an energy theorem. The new state is point based and requires
an explicit compatible cylindrical/field representation before adoption.

## Continued common-radial-face controls

Seven additional remote 3f4120fa cases completed, passed original-platform
bitwise/physical audits and were retrieved to `outputs/face_radial_remote_controls`.
Together with the five local runs this gives **12 of 14 completed/audited cases**;
the nr16/32 full contacts remain assigned to the live checkpointed worker.
Consult its saved process observation before considering any restart.

The nr32 acoustic density/electron-temperature errors are
1.69024e-8/1.12770e-8, continuing the nr8/16 refinement. Its physical energy
drift is −8.76383e-12, with finite-step remainder −5.83e-14. The nr8 quarter-step
contact retains −4.10007e-9 physical energy drift while its finite-step remainder
drops to −2.53e-15, confirming that prototype's spatial transport limitation.
Direction/quadrature/time controls and nr32 rotation are also preserved.
Maximum independent local-energy agreement in those seven audits is 3.13e-15.

## Next action and retained scope

Extend the derived joint flux to the actual cylindrical physical state, with
regular-axis forces, angular momentum and true radial kinetic/thermal metrics;
test original full two-dimensional equations and unforced dynamics. Then prove
and verify the shared density/current/magnetic energy and species enthalpy/
Poynting coupling. The existing physical vacuum, heated edge, transport/
ionization and formation requirements remain unchanged. Both prototypes remain
verification components; neither proves FRC formation, confinement or stability.

The inspected figure is `outputs/two_entropy_diagnostics/report.png`. Sources,
raw arrays, complete accepted histories, quarter-run checkpoints, physical
audits and separated species diagnostics accompany the reported results.
