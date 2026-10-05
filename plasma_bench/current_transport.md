# Physical current advection and drift-stress candidate

Uses derivatives of actual reconstructed fields and both one-sided traces.
The volume advection form is skew in current/test-current plus the actual
changing-density half term. Central skew face terms retain normal-current
jumps. Radial wall/axis numerical transport flux is zero; periodic axial faces
pair. Actual reconstructed boundary traces are not claimed identically zero.

Bulk drift stress is the velocity-gradient adjoint with shared central radial
tractions. The electric stretching term includes the exact paired velocity-jump
work. Symmetric stress preserves total axial force and actual angular torque.
Self-drift advection uses the same skew form with w=J/(en), zero carrier-density
rate, and the exact physical mass fraction. Its consistency still requires
measuring the reconstructed strong/current-jump continuity defect.

`pair` returns separate advection, stretching, self-drift and bulk stress terms.
It adds no ideal/Hall, species pressure, thermal source or time evolution.
Physical equation consistency, near-axis/wall rates and interface quadrature
must be verified separately from the algebraic work identities. This candidate
must not replace a verified source until those gates pass. No diffusion or
energy repair is introduced.

## Verified component scope (2026-10-04)

The unchanged candidate now has sixteen frozen cases, bitwise replay and an
independent continuous-field audit. See docs/CURRENT_TRANSPORT_VERIFICATION.md
for convergence, local/axis/wall errors, charge defects and the measured grid/
continuous energy difference. RegularPlasma now combines this unchanged component
with both species and magnetic evolution; full fast-wave periods have independent
verification. Long contact, local interface budgets and formation remain open;
see REGULAR_PLASMA_VERIFICATION. No production formation solver was replaced.
