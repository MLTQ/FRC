# Verified implicit wave evolution; unresolved nonuniform-contact growth

The full regression suite passes **432 tests in 227.47 s**. Seven new focused
checks cover independent uniform/moving linearization, real collocation
equations, original strict residuals, actual returned stages, preserved inputs
and atomic exhaustion. All simulation, replay, audit, report and test jobs in
this checkpoint are terminal. Hardware remains prescribed waveforms; no remote
or publication was used. These are filled periodic verification cylinders, not
the user's proposed quartz-tube formation experiment.

## Solving the original time equations

Physical block preconditioning, implicit midpoint, fourth-order Gauss and
third-order Radau evolve the complete actual fluid/current/species/magnetic
state. Equations, scales, iteration, numerical temporal damping and limitations
are in IMPLICIT_PLASMA_NUMERICS.md and module companions. No original spatial
term, physical mass or induction residual is weakened or removed.

Initial chord-only midpoint/Gauss pilots rejected their long contact before a
complete transit. The saved Gauss eighth accepted state reproduces a ninth-step
chord stagnation at 9.105566e-12, above a 5.403121e-12 nonlinear limit. A
matrix-free Newton correction reduces it to 3.027969e-14 under the same bound,
with ten Krylov iterations. That captured input is a regression test. The old
source archives and failures remain unchanged; successful nonlinear iteration
does not establish physical trajectory accuracy.

## Eleven complete original-parent waves

All initial/final/penultimate arrays, step residuals/limits/iteration counts and
non-runtime metrics replay **bitwise for all eleven complete trajectories**.
Each has independent initial/final physical quadrature orders 16 and 20, full E,
both temperatures, B/J and separately reconstructed continuous energy. Physical
argon/electron masses, amplitude 1 m/s in equal-energy coordinates, two degree-5
radial fluid elements, 16 axial points and the existing field grid are fixed.
The period is approximately 116.5 ns. Nonlinear time tolerances are explicitly
atol=1e-15, rtol=1e-12 for this small wave, with original induction tolerance
1e-10. These are not selected laboratory parameters.

| Method | Complete steps/period | Successive physical B difference ratios | Time order |
| --- | --- | --- | --- |
| Midpoint | 30, 60, 120 | 3.975482 | 2 |
| Gauss | 30, 60, 120, 240 | 15.967079, 15.991775 | 4 |
| Radau IIA | 30, 60, 120, 240 | 7.975201, 7.994196 | 3 |

These ratios use independently reconstructed physical B differences, not merely
coefficient norms. The total continuum B error approaches about 2.224e-3 on
this deliberately coarse fixed spatial mesh; time refinement does not remove
spatial/nonlinear error. Gauss/midpoint continuous energy changes divided by
wave energy fluctuate around 1e-7 or less at this weak-signal precision level;
they do not show exact conservation or clean time convergence of that tiny
energy diagnostic. Grid and continuous values remain distinct.

Radau's measured physical energy loss divided by initial wave energy is
1.592348e-3, 1.999813e-4, 2.505569e-5 and 3.049647e-6 at 30/60/120/240 steps.
This is explicit numerical temporal damping, not a physical transport closure.
Its convergence is measured without any energy/temperature correction.

Sources: Gauss/midpoint trajectory closure
`feeb358468dfd0074e61effcca9af8dcb727748f113907d22b4e0e74e958c8ef`;
Radau trajectory/audit closure
`445623fdadd9f0825ae5441dc4c3fe26e92859e5c5d53564d5639949c9bb4240`.
Raw records are outputs/implicit_wave_{method}_{steps} and
outputs/radau_wave_{steps}, with separate *_audit records and source archives.
The inspected PNG/PDF report and exact measured table are in
outputs/implicit_plasma_report.

## Long contacts did not validate

Three Gauss controls and a Radau control were stopped after deviations from the
exact contact grew substantially. Every accepted nonlinear/induction equation
still passed its original bound. These are **incomplete and inaccurate runs**,
not completed transits or reported nonlinear-solver failures. Their termination
records preserve the last saved accepted checkpoint; it may precede additional
accepted states in memory at interruption. No unsaved state is claimed audited.

| Run | Last saved time | Independent density signal error | Independent max unwanted B |
| --- | ---: | ---: | ---: |
| Gauss, 4 elements / 120 target steps | 50 microseconds | 0.551685 | 1.341993e-3 T |
| Gauss, 4 elements / 240 target steps | 30 microseconds | 1.351799e-3 | 2.697311e-5 T |
| Gauss, 8 elements / 120 target steps | 20 microseconds | 2.902119e-2 | 6.069165e-4 T |

The three Gauss checkpoints have separate order-16/20 field, temperature and
continuous-energy audits in outputs/implicit_contact_failure_audit. Radau's
last saved 30-microsecond state also has growing magnetic energy, so L stability
alone does not resolve this problem. It is preserved in
outputs/radau_contact_ne4_t120 and is not claimed to have the same independent
checkpoint audit. No run completed an accurate requested 100-microsecond transit.

## Spatial growth exists without a time integrator

Setting the contact's uniform axial velocity to zero defines a separate exact
stationary continuum state. Three complete 360-coordinate centered-difference
spatial matrices are frozen and replay bitwise in
outputs/stationary_contact_linear_probe. The uniform rest control has maximum
real eigenvalue 4.28e-9 /s; the nonuniform contact has 854532.845 /s and
854532.963 /s at perturbation amplitudes 1e-4 and 1e-5. Corresponding nonzero
mixed-direction errors decrease from 5.637437e-8 to 6.137574e-10. All original
strict induction bounds pass. Source:
`b82a7bf386a9d591e67699db71c9a5326e35c4bf709a30364b208645edd2449f`.

The uniform leading vector is a numerical null mode: its original relative
mixed-direction error is ill-conditioned and is explicitly retained. A proper
nonzero independent rest-operator control instead gives errors at most
4.402961e-12. Do not present the null-mode relative value as a useful accuracy
test. Component ablations reduce the 854533 /s rapid growth to at most about
0.153 /s when the paired species-pressure/current-entropy feedback is removed;
fixed-bulk species/field and magnetic-only blocks have maxima near 1e-6 /s.
These incomplete diagnostic matrices are not physical evolution models or a
justification for dropping the species pair.

## Separate continuum linearization and remaining interpretation

An independent original-equation Fourier/Bessel linearization imports none of
the production spatial operators. Its three uniform finite-mass wave branches
have equation residuals 7.26e-13, 5.41e-15 and 7.03e-15. Eight variable-profile
matrices/leading vectors and all uniform metrics replay bitwise. At the first
Bessel radial mode, largest real growth converges from 13061.78 /s (8 axial
points), to 8759.13 (16), to 8763.749276 (32), to 8763.749254 (64/128).
Modes 4/8/16 at 64 axial points give 11815.41/12157.90/12307.16 /s. Leading
magnetic-divergence relative residuals are <=4.48e-10; eigensolve backward errors
are <=8.55e-17. Source/audit/replay are in outputs/stationary_contact_audit,
source `b7dde0bbad1b307e752bc2a9b441e158dd7c790dd88728c684782667aee647a2`.

This independent calculation also has growing modes. Its radial spaces differ
from the production matrix; the rates are **not a matched-mode comparison**.
Do not conclude that all contact growth is numerical, or that the rapid
production mode is physical. The Euler density-wave stability literature
([Ranocha and Gassner](https://arxiv.org/abs/2009.13139)) motivates this diagnostic
but is not evidence about this plasma's particular instability.

Next, under frc-7mm.7: reconstruct the production leading mode in physical
fields, compare matched continuum radial/axial subspaces and independent
perturbation evolution, and measure refinement. Resolve any demonstrated
coupling defect and the original local enthalpy/Poynting interface budget
without fitted states or energy repair. Then repeat complete long contacts.
Only afterward proceed to heated moving exact vacuum, physical transport and
ionization, coil/end boundaries, converged formation and the experiment UI.
The full FRC goal remains open; no 3-D stability or confinement is established.
