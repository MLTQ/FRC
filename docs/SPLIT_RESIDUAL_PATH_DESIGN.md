# Split-inverse residual paths and component sensitivity

This continuation follows the completed checkpoint 5ed1d8e. That turn made
verified progress; all its workers are terminal. The overall FRC goal stays
open. Verification uses complete saved-input workflows, independent physical
reconstruction and fresh full replays, with no unit tests.

## Exact rejected-state observation

Use the preserved F64/256, B256/256 minus input at epsilon 1e-4 and tolerance
1e-11 in radial16_finest_controls/case_01.npz. Reconstruct and require bitwise
identity of the base state, perturbation and actual signed state/field. Assemble
the original full fluid, ideal/Hall, current-transport and species forcing with
the physical masses and actual density. Use the identical direct-metric setup.

Let T be the admissible lift, D = H T its stored dual lift, F the original
wall curl, A the electric mass inverse and I the physical current inertia.
The original split-inverse core evaluates

    b = T y, j = -A F^T H b, e0 = A f, ei = A I j, e = e0 + ei,
    rhs = D^T F e0,
    K y = T^T D y + D^T F A I A F^T D y,
    u = K y - rhs, v = D^T (F e - b).

In exact arithmetic u + v = 0. The observed finite-precision gap may include
H/T association, current and inertia actions, electric addition, curl and
dual-projection distributivity, and outer vector arithmetic. Save the vectors
and each of those contributions; retain any unexplained remainder. A large
term norm alone is not proof that it explains the vector discrepancy.

An observation-only copy adds records after evaluating the original residuals.
It keeps every original operation in CG/refinement/acceptance, strict bound,
iteration budget and failure message. Its first-attempt rejection must match
the archived original exactly. The complete five-attempt original rejection
already has a separate full replay; the observation does not accept a full
plasma rate. Repeat the complete reconstruction and observation independently.

## Physical component audit

Reconstruct the accepted signed rates at epsilon 1e-4 and 1e-5 at all original
quadrature 16/20 points. Require the original physical fingerprints before
interpreting component error or sensitivity. Use the independent radial16,
N128 mode, component reference norms and fixed physical axis/wall bands.
Retain radial profiles and whole-field norms rather than inferring local
accuracy from the dominant toroidal-current component.

## Unchanged physical scope

These are filled, periodic argon verification fixtures with an impermeable
radial wall and prescribed constant bias. Parent equations and assumptions
remain those in RADIAL_MODE_DESIGN.md. No numerical diffusion, filter, density
floor, artificial mass, atmosphere, heat or energy repair is introduced.
Spatial/truncation and floating-point errors are measured, not reclassified
as physical transport. Heated moving exact vacuum, physical transport and
ionization, real ends/coils, complete trajectory convergence, formation and
the experimental UI remain required. The user's lab envelope is provisional.
