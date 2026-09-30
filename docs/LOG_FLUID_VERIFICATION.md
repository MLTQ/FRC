# Self-consistent fluid diagnostic checkpoint

This is a progress checkpoint, not adoption of a full plasma solver. The complete
local suite passes **376 tests in 151.46 s**. All nine two-dimensional parent-rate
cases and eight local trajectory records (including two preserved failures)
have independent physical audits. The 2720×1600 six-panel diagnostic figure was
visually inspected, including corrected resolution tick labels.

The original fourteen-case remote baseline has finished: nine complete cases
and five failed contacts, all fourteen records audited on their original platform
and retrieved with raw-array hashes verified. The separate fourteen-case remote
radial-work comparison is still running. Its incomplete results are not treated
as final evidence here. Exact process handles and next actions are in GOAL_STATUS.

## Sources, model and tests

Baseline source: `0aa417e588be92ab1ad5d1d5e1deeece7e6c2f02cb54ded16a01562849828fea`.
Radial-work source: `4f0c824a662e581dc78ddfbd81a2c9cc06342fb85d19d62fda44d0fb709ab3a9`.
Both archives retain dependency file hashes, source copies and library versions.
Previous bulk/current/thermal and older nonlinear calculations are unchanged.

The model has actual n, mr/mz/angular momentum, two physical log entropies and a
separate conservative energy account. Domain: R=.035 m, z=±.12 m, periodic z and
impermeable regular-axis/radial walls. Physical masses include electron plus
39.948-amu ion mass; n0=1e18 m^-3 and total reference pressure .4 Pa are verification
inputs. There are no magnetic/current dynamics or external material. No
numerical diffusion, energy reset, density/pressure floor or hidden retry is used.
Equations and work derivations are in LOG_FLUID_NUMERICS/RADIAL_WORK_NUMERICS.

Independent continuum checks include finite-amplitude characteristic waves,
opposite propagation, a translating two-temperature contact, rotating equilibrium
and fully 2-D Euler/pressure rates. Physical velocity/curl and exact K derivatives
are checked separately. Atomic failed stages leave the accepted state intact.

## Preserved failures and their diagnoses

The original nr8 contact fails at attempted step125 of516, after0.24031 transit,
with nonpositive reconstructed face density after artificial radial motion.
Its accepted state/history/traceback remain in `outputs/log_fluid_pilot`.
The focused test suite did not establish long-duration stability.

Seven independent acoustic spectrum controls show degree5 nr8/16/32 growth rates
5.87362e5, 7.45952e5 and1.34395e6 s^-1. Original nonlinear RHS central differences
verify the linearization. An independent matrix exponential and RK4 time
refinement retain the growth; it is a spatial defect. The energy-adjoint linear
control has maximum real part below6e-10 s^-1 in those same cases.

Nine translation probes show a separate obstruction: for m_z=U n_cell, any
globally conservative force has zero kinetic power, but the current thermal
update has a convergent chain-rule defect. Electron/ion thermal powers atnr8
are5.08993e-4/4.98473e-4 W, dropping to2.10329e-9/2.10367e-9 W atnr64.
Reverse flow reverses the sign and higher quadrature agrees. A force-only exact
energy correction would therefore violate momentum conservation. No such
correction was applied.

The first radial pressure-only adjoint trial spoiled centrifugal balance.
Its source/measurements are preserved. Pairing both pressure and kinetic
convective/centrifugal work reduces the rotating radial-force RMS to
.328377/.0126854/.000423483 m/s² atnr8/16/32. This is not sufficient for adoption:
the fully 2-D force comparison below exposes a separate local defect.

The remote baseline also demonstrates that simply completing a wave is not
convergence: density relative L2 errors atnr8/16/32 are4.84313e-5,1.00744e-6,
3.29836e-5. Finer radial resolution eventually amplifies the unstable modes.
All five remote contacts reject, including reversed flow and half time step.
The exact rejection step differs slightly by platform; original platform/library
provenance and replay are retained rather than claiming cross-platform identity.

All five actual remote rejected steps replay their original error and leave the
retained accepted arrays unchanged. Their already unstable last states expose
an additional quadrature limitation: order12 versus order22 local physical-energy
differences reach1.05e-3 relative, and order18 versus22 still differs by1.73e-4
in the worst failed case. Those energy integrals are not converged. The original
audit accuracy assertion and log were preserved; the expanded audit reports this
warning separately from successful state/history/rejection replay. No solver
tolerance or accepted result was changed. Both local failed-state checks remain
within their stated quadrature tolerance.

The revised nr8 contact initially fails at step30/516. Its higher radial acoustic
frequency makes the original coarse fixed step marginal at variable sound speed.
Explicit half- and quarter-step runs both complete one full transit; their final
physical density/velocity/Te differences are6.29e-11/9.55e-11/8.79e-11 relative L2.
The original larger-step failure is retained in `outputs/radial_work_pilot`.

## Completed local contact and physical energy

For the revised nr8 half-step full transit:

| Quantity | Measured value |
|---|---:|
| Density relative physical L2 error | 2.90525e-5 |
| Electron temperature relative physical L2 error | 5.15852e-5 |
| Ion temperature relative physical L2 error | 5.29554e-5 |
| Physical energy relative change | −4.10015e-9 |
| Conservative energy-account relative change | −1.45e-16 |
| Integrated spatial energy-defect / initial energy | −4.10007e-9 |
| Remaining finite-step energy defect / initial energy | −7.36e-14 |

At quarter step the remaining finite-step defect is−2.67e-15. Physical energy
is not mislabeled exactly conserved, and the account is never used to repair it.
Particle and axial-momentum changes in both controls remain within5.56e-16
relative. Initial projection errors and all five saved times are separate from
final continuum error in the raw records.

## Remaining first-annulus radial-force error

Nine original 2-D parent-rate cases cover degree3/5, nr8/16/32/64 and increased
quadrature. All complete state/rate arrays replay bit for bit. Independent
order18 physical local-energy rates agree within5.24e-15 relative.

| Degree5 radial cells | Full acceleration relative L2 | First-annulus radial component relative L2 |
|---:|---:|---:|
| 8 | .501467 | 4.95635 |
| 16 | .109340 | 6.33007 |
| 32 | .0285246 | 6.65104 |
| 64 | .00734626 | 6.73205 |

The first-annulus axial/azimuthal component errors atnr64 are8.54e-10/3.97e-12.
Thus the local radial error cannot be excused by shrinking axis-cell volume or
small global energy balance. Near-axis compatibility among actual momentum
means, the divergence-commuting particle lift and kinetic/convective work is
the next unresolved derivation. The cause is not yet fully established.

A four-mesh force-term diagnosis further isolates the problem. Atnr64 the
first-annulus pressure-source and kinetic/centrifugal-source relative errors are
6.44600 and6.44732; conservative convection is2.15e-10. Large compensating source
errors can preserve rotating equilibrium, so equilibrium residual alone is not
a force-accuracy test. Raw individual forces and analytic targets are retained
in outputs/radial_force_terms. The next correction must address these two work
adjoints without spoiling the already accurate convective transport.

## Audit and continuation

All eight local trajectory records pass initial-state/rate replay; the six
completed trajectories also pass last-step bitwise replay. Failed trajectories
remain failures. Independent order18 physical local energy at every saved
snapshot agrees within1.04e-15 relative; recorded history totals replay.
Raw sources, fields, every accepted global history, ledgers and tracebacks are
under `outputs/log_fluid_*` and `outputs/radial_work_*`.

Finish the original-platform remote audits when each fixed-source worker exits;
do not restart a live case or alter its closure. Then resolve the local radial
force defect before replacing the full moving state. Magnetic/current coupling,
physical free edge/exact vacuum, transport/ionization, actual end/coil boundaries,
shocks and converged formation remain open. Axisymmetry does not prove 3-D
stability or particle confinement. Work remains local/private and unpublished.
