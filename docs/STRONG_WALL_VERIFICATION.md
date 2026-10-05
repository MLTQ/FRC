# Bounded-wall correction and verification

The full experiment sandbox remains unfinished. This checkpoint identifies a
large wall-space defect, implements an explicitly boundary-conforming candidate
and verifies complete unforced contact transits and fast waves. **436 tests pass
in 232.03 s.** These are filled, periodic, fully ionized verification fixtures;
the wall is not a model of the proposed quartz tube, vacuum edge or open ends.
Read STRONG_WALL_DESIGN.md for equations and boundary assumptions.
All new jobs are terminal. Every successful trajectory and diagnostic case
has full bitwise replay. The visually inspected scientific report is
[convergence.png](../outputs/strong_wall_report/convergence.png) and
[PDF](../outputs/strong_wall_report/convergence.pdf); its JSON retains exact
time ratios, both quadrature orders and all source/input identities. All 32 new
source archives verify, 22 earlier calculation closures are unchanged, and
32 archives from the preceding checkpoint are byte-for-byte unchanged.

## Identified rapid mode

Independent root-polynomial reconstruction of the frozen 854532.962678 /s mode
finds outer-wall radial velocity and Btheta RMS about 5.36944 and 4.99716 times
their respective volume RMS. The outermost radial field cell contains 95.18%
of radial-velocity squared norm and 79.10% of toroidal-B squared norm. Normal
wall current RMS is 0.254699 of total volume-current RMS. These physical traces
violate the originally intended impermeable, Btheta(R)=0 boundary despite
passing the old weak equation checks.

The strong reconstructed Ampere discrepancy is 1.42923 relative to current;
strong div J / current-gradient norm is 0.0571344. These are not complete
distributional interface budgets. The coarse interpolated background also has
ion-pressure error 0.00441230 and total-pressure maximum error 0.00159979.
Fourier/Bessel projections explicitly measure uncaptured norm. At 32 radial
modes, transverse flow and Btheta still have substantial uncaptured norm because
the original wall traces are incompatible with that continuum basis. The old
mode is not declared a perfectly matched continuum eigenfunction.

Applying separately derived physical K-metric Galerkin restrictions to the
saved original matrix gives largest real rates:

| Restriction | Growth, 1/s |
| --- | ---: |
| Original matrix | 854532.962678 |
| Normal fluid velocity only | 335136.225099 |
| Toroidal magnetic trace only | 14601.194205 |
| Both intended wall traces | 9203.373408 |

The actual full constrained RHS is differentiated independently in 344 active
coordinates. Contact difference amplitudes 1e-4 and 1e-5 give 9203.372137 and
9203.373408 /s; mixed-direction errors are 5.47878e-6 and 5.22837e-8. Its matrix
agrees with the separately restricted original matrix to 1.12125e-11 at the
smaller amplitude. Uniform rest has maximum real part 3.00507e-6 /s on an
oscillatory mode at 5.32563e6 /s, with independent mixed-direction discrepancy
2.48605e-10. It is not evidence for physical uniform-state growth.

All three full matrices and every physical/projection array replay bitwise.
Maximum strict induction residual is below 3e-11; normalized toroidal-rate
wall traces are below 1.2e-16. An earlier larger-space iterative solve rejected
a small trace leak and is preserved in outputs/contact_wall_diagnosis with a
terminal failure record. The accepted implementation uses an admissible Krylov
basis, retaining full unreduced checks; it does not relax the trace tolerance.

Source: ded3bac43979f5bde181c2eb6989d6af28b57f1d64bca09dd3f94413b7203a9a.
Evidence: outputs/contact_wall_diagnosis_conforming, its independently generated
replay, and outputs/contact_wall_replay_verification.

## What changed mathematically

The new declared variational curl is Fw=P F, where P is magnetic-energy
orthogonal projection onto the exact toroidal wall-trace null space. All time
stages and inner Krylov vectors are parameterized in this space. Normal fluid
wall momentum is inactive; tangential velocities remain free. Physical masses,
both evolving species temperatures and actual density-dependent inertia remain.

Full weak Ohm, Ampere and Fw Faraday equations retain strict checks. The old
unprojected F E-Bt is reported separately and is **not** called zero or a passing
old-incidence equation. Its magnetic-energy dual vanishes on admitted tests.
Independent physical trace, reaction-work and actual energy-derivative tests
verify the derivation. No atmosphere, floor, viscosity, time damping, state
filter, fitted evolution field, altered mass or heat/energy repair is added.

The cell-integrated analytic toroidal initial field has a declared one-time
boundary-conforming approximation. Its error is measured separately. The
Gauss integrator uses the same fourth-order full nonlinear stage equations;
its fixed-background restricted inverse is only a solve preconditioner.

## Complete trajectories and conservation

The final source for the conforming trajectory protocol is
d7a422b0f3e147e84e1acbb5b3db5c17b4eb162eab5d2d55845d829e18ebbc96.
Complete-array and non-runtime-metric replay and independent physical
quadrature16/20 are required for every reported run. Full E, both temperatures,
point/sample positivity and axis/wall/max field errors remain explicit.

Wave time differences at 30/60/120/240 steps give physical B ratios 15.967079
and 15.991774. Contact density differences at 60/120/240 steps give 15.883252.
These are differences between solutions at fixed spatial resolution, separate
from total analytic-reference error. Contact physical relative energy drifts
at those steps are -3.044159e-8, -1.905815e-9 and -1.191590e-10, consistent with
fourth-order temporal convergence. Grid and independent physical energy agree
to the measured precision. Integrated conservation is not a local face budget.

At eight radial fluid elements, the full fast-wave B error is 7.31894e-5,
full-E error 2.53448e-5 and density/pressure error about 0.00561855. The
two- and four-element cases use the same eight-cell field grid, retaining the
same magnetic error floor; this is recorded, not claimed as mesh convergence.
The wave energy change relative to its perturbation energy is 4.63884e-8 at
eight elements. Tiny wave-energy drifts retain a floating-point precision floor.
The raw-incidence curl defect at that final wave is 3.75696e-6 relative to the
magnetic rate in the reported coefficient norm; it is not an acceptance residual.
Near-zero contact B rates make this relative defect ill-conditioned, so the
absolute unwanted B/current and independent original E must also be read.

The pre-basis-change four-element contact also completed and replayed bitwise,
with source 20a47715cfbb665541a05fa332ba139b989eceacdce17d40244a07b110c9c275.
Its original source is in its archive; the updated source is deliberately a
different identity. This successful pilot does not erase the separate rejected
mode probe.

The three complete 120-step contact transits give these order-20 physical errors:

| Fluid elements / axial points | Density signal error | Full E error |
| --- | ---: | ---: |
| 2 / 16 | 0.00201835 | 0.0210529 |
| 4 / 32 | 0.0000506225 | 0.0104219 |
| 8 / 64 | 0.000000997995 | 0.000500415 |

Finest electron/ion-pressure errors are 3.77238e-6 / 3.86838e-6. Unwanted B is
1.00896e-13 T and current 1.12468e-4 A/m² maximum at independent order20; neither
is silently set to zero. Particle relative drift is 1.49747e-17, and entropy
changes per initial particle are -1.71456e-17 / -1.47257e-18. Point and sampled
densities and both physical temperatures remain positive. The finest contact
uses three Newton corrections under its unchanged acceptance bounds. These
small reference errors do not establish the separate perturbed-mode accuracy.

## Independent growing continuum mode

The frozen original-equation Fourier/Bessel N=64 first radial mode grows at
8763.749254 /s. Its full electric field is derived from ideal, Hall, physical
inertia and both pressure gradients, independently of production projection.
Original Faraday holds to 1.42219e-9 relative; its small magnetic divergence
residual is retained. The analytic J1 wall conditions hold to floating precision.

Sample this same continuum mode into the actual conservative state and analytic
magnetic flux/cell toroidal averages, then evaluate both signs of a small
perturbation through the full nonlinear RHS. Centered *physical* rates, including
the differentiated velocity/entropy quotients, are compared with lambda times
the original fields. No production eigenmode is fitted. All twelve mesh,
amplitude and quadrature cases and metrics replay bitwise.

| Fluid elements / field cells per direction | Full E error | Radial velocity rate | Toroidal B rate | Current rate |
| --- | ---: | ---: | ---: | ---: |
| 4 / 8 | 0.286066 | 0.590992 | 32.0863 | 27.6674 |
| 8 / 16 | 0.0833377 | 0.051017 | 4.08666 | 5.56066 |
| 16 / 32 | 0.00648538 | 0.002272 | 0.256271 | 0.480471 |

The finest radial-B rate error is still 0.547444. Initial magnetic representation
errors are much smaller, and slow growth comes from cancellation of faster
terms: accurately representing B alone does not certify B_t. The epsilon1e-4
versus1e-5 check at eight elements barely changes these errors, supporting a
spatial rather than finite-difference-amplitude limit there. Quadrature16/20
agrees; **these mode-rate errors do not yet meet a formation adoption gate**.

Source: c100c3eea05acbb767c65887bf692335ff88b788d884eef3839aaf7e5b0ec4fb.
Evidence: outputs/continuum_wall_probe, its independent replay and
outputs/continuum_wall_replay_verification. These are instantaneous directional
rates, not completed independent perturbation trajectories or laboratory
stability/confinement predictions.

## Required next work

1. Refine radial and axial field resolution separately on the independent
   nonuniform growing mode; retain full physical rates, current and E, initial
   errors, axis/wall/max norms and quadrature controls. Resolve the slow magnetic
   response rather than only checking the faster-wave temporal method.
2. Evolve small continuum perturbations independently with amplitude/time/mesh
   controls. Preserve all boundaries and physical masses; compare actual fields
   and growth, not just the largest eigenvalue of differing spaces.
3. Close local physical charge, enthalpy and Poynting interface balances,
   including distributional jumps. Retain original cold-flow/weak-pressure
   limitations and avoid promoting integrated energy cancellation to local proof.
4. Then proceed through heated moving exact vacuum, transport/ionization,
   actual coil/end conditions, converged formation and experimental controls.
   The whole frc-goal and current numerical acceptance beads remain open.

The user's provisional one-metre tube, 10–20 cm diameter, selectable gas,
0.1–100 mTorr neutral fill and adjustable Helmholtz-like bias remain exploration
inputs, not selected operating recommendations. Hardware is prescribed waveforms.
