# Common radial momentum verification checkpoint

The complete suite passes **382 tests in 169.23 s**. Six new tests cover an
independent regular polynomial force and physical kinetic energy, the retained
diagnostic work expression, original two-dimensional acceleration/curl refinement,
unforced nonlinear waves, atomic failure, actual pressure/thermal work, and
bitwise recovery of an interrupted trajectory. The new component remains
experimental; see FACE_RADIAL_NUMERICS for equations and energy limitations.

Source identity:
`3f4120faf1b96ee2e94453339efe215399071fd632827aa6ae459cd9c5126737`.
Old fluid 0aa417e5, radial-work 4f0c824a and earlier verified source closures
remain unchanged. Full hashes, package versions, inputs, arrays and source
archives accompany each study.

## Independent original-equation rates

All nine degree3/5 mesh and quadrature cases are reproduced bitwise and audited
with independent order18 physical quadrature. Local physical energy-rate
agreement is within 5.25e-15. The degree5 results are:

| Radial cells | Full acceleration relative L2 | First-annulus radial relative L2 | Absolute scaled global energy-rate defect |
|---:|---:|---:|---:|
| 8 | 1.98671e-4 | 1.75954e-5 | 8.28819e-7 |
| 16 | 2.69154e-6 | 1.91939e-7 | 1.40178e-8 |
| 32 | 5.60506e-8 | 3.08493e-9 | 2.23570e-10 |
| 64 | 1.32683e-9 | 1.19443e-10 | 3.51141e-12 |

The energy-rate scale is initial physical energy times sound speed/radius.
Its nonzero value is a retained spatial transport defect, not solver tolerance
or a repaired temperature. The older radial-work method's first-annulus radial
errors are 4.95635, 6.33007, 6.65104 and 6.73205 for the same physical state.

At nr64 the new independent first-annulus pressure, centrifugal and advection
relative errors are 1.23224e-10, 5.41147e-13 and 1.34151e-10. Individual whole-
domain pressure/centrifugal response errors remain 0.01991/0.08343, because
each individual force violates the imposed outer-wall constraint although their
sum satisfies it. They are retained in the records, not confused with the
old inconsistent near-axis source pair. The independent wall-compatible
polynomial force is reproduced at approximately floating-point accuracy.

The preliminary common-face prototype with an exact radial total-work correction
had full acceleration errors .00243039, .000612014 and .000190902 at nr8/16/32.
Its failed high-order convergence check is preserved. The adopted experimental
derivative uses the physical centrifugal force and records convective energy
error rather than changing that force to cancel it.

## Stationary radial acoustic limit

Seven degree/mesh/quadrature controls retain complete matrices/eigenvalues.
Central differences of the actual nonlinear RHS agree with the proposed linear
operator at three amplitudes. Its physical energy adjoint defect is near
roundoff (at most 3.07e-16); it has no resolved exponentially growing acoustic
eigenmode (maximum real eigenvalue times R/c is 5.52e-15). At perturbation
amplitude 1e-5, the largest original-RHS difference error is 8.04e-10.

Independent cylindrical frequencies use the first three positive J1 zeros
with regular J0 pressure and zero wall-normal velocity. At degree5 nr32 the
absolute relative errors are 1.01e-15, 1.67e-13 and 4.62e-12. The fundamental
frequency has reached floating-point resolution; this is not evidence of
comparable accuracy for all represented short wavelengths. Independent matrix
exponentials preserve two-period physical linear energy near roundoff, and
RK4 time halving reduces the selected mode error by approximately sixteen.
These checks establish the stationary acoustic limit, not nonlinear stability.

## Completed local trajectories

Five unforced cases finish and pass initial/rate/last-step bitwise replay,
independent order18 physical-energy integration and complete saved-history
checks. Maximum local energy agreement error is 1.37e-15. All five remain
positive without a floor or state repair.

The nr8 contact finishes one full 195.629-microsecond transit using the declared
half-sized step. Final density/electron-temperature relative errors are
2.90525e-5/5.15852e-5. Physical energy drift is −4.10015e-9, the conservative
account changes by 4.35e-16, and the separated finite-step defect is −7.41e-14.
Most of the physical energy drift is therefore spatial, not time-integration
error. This agrees with the earlier translation obstruction; fixing radial
force accuracy does not establish exact axial nonlinear energy transport.

Acoustic nr8/16 density errors are 4.84313e-5/1.00740e-6; electron-temperature
errors are 3.21717e-5/6.71694e-7. Physical energy drifts are −4.52186e-8 and
−5.99532e-10; remaining finite-step energy defects are −1.90e-11/−1.02e-12.
Rotating nr8/16 density errors are 1.27e-10/2.30e-12; their energy changes are
at floating-point resolution. Initial projection, maximum norms and physical
small-signal rate measures remain separately available in every record.

The remaining finer-grid and time/sign/quadrature controls are assigned to the
checkpointed study on the authorized remote host. Do not count a running
checkpoint as a completed or audited trajectory. See GOAL_STATUS for current
worker state. The scientific six-panel report reads only completed local
records: `outputs/face_radial_diagnostics/report.png`.

## Recovery of the previous remote study

The old 4f0c824a worker stopped after its SSH connection reset; a subsequent
process inspection found no Python worker in that study directory. Thirteen
saved records were audited on their original platform and retrieved:
**12 complete, one numerical failure**. The failed attempted step replays with
unchanged accepted arrays and no higher-quadrature warning. The nr32 contact
has no saved terminal record and is explicitly **interrupted**, last logged at
75% after 7198.1 s. No final state can be reconstructed from that log.

These records are in `outputs/radial_work_remote_recovered`, including the
source archive, hashes, original-platform audits and recovery status. The
missing old case is not rerun because its force formulation is superseded and
known to be inaccurate; it is never counted as complete. The new checkpoint
driver addresses the lost-computation problem without modifying old records.

## Next requirements

Finish and independently audit the remaining unforced controls. Derive shared
axial/radial convective momentum and kinetic-energy transport that satisfies
the original equations and physical conservation simultaneously. Preserve the
momentum/translation obstruction rather than applying a force-only energy fix.
Then integrate magnetic/current/drift stresses and complete species enthalpy/
Poynting work, followed by the original heated-edge/exact-vacuum requirements,
physical transport/ionization and converged formation. Neither these bulk tests
nor an axisymmetric closed-flux fixture establish 3-D stability or confinement.
