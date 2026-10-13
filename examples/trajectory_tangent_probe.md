# Actual spatial tangent along an observed trajectory error

`run(destination, trajectory, verification)` requires a terminal refined
F8/64,B32/64 trajectory and its complete byte-identical state/stage/failure
replay. It re-creates the original split-inverse fixture and checks initial
fluid/magnetic bytes against that exact initialization. Step10's complete
conservative error subtracts the declared independently growing linear tangent
from the saved change, using expm1 and all components. No mode is discarded.

Both this saved error and the independent original tangent are pulled back to
the existing admissible ModalWallCoordinates and normalized by their coordinate
infinity norms. Native errors, fixed scales, radial null basis, lifted tangents,
pack/unpack differences and constraints are archived. This is a diagnostic
input construction; no accepted state is projected, filtered or overwritten.

The complete original stationary spatial RHS is evaluated at both signs of
each normalized direction, amplitudes1e-4,1e-5,1e-6. The physical masses, both
species, current transport, wall space and strict original1e-10 induction
checks remain. This module constructs no time integrator or dense Jacobian.
The discrete base rate, including its stationarity defect and actual input
immutability, is measured first. Each centered magnetic action also passes
through the unchanged production weak Ampere map for a separate current-rate
consistency comparison in the physical audit.

Every sign's inputs are archived before the actual call. Their in-memory
state/field fingerprints are checked afterward, including a caught RHS
exception. Accepted full fluid/magnetic rates, current rate, solved E,
positivity and all diagnostics are saved immediately. A returned exception
is recorded verbatim, preserves any preceding sign and continues to the next
declared amplitude; it is not converted to a rate or automatically called a
physical/numerical instability. Other workflow errors propagate with records.

Centered actual conservative rates give the original coordinate action.
Projection onto the input and its orthogonal remainder use the declared fixed
coordinate metric; the projection is not an eigenvalue or physical growth
claim. Separate independent quadrature16/20 physical reconstruction, continuum
mode comparison, amplitude differences and complete fresh replay are required.
Physical-RHS differentiation and the conservative-Jacobian pushforward must
remain distinguished when the discrete base rate is nonzero.

Run through the guarded tangent workflow with8GiB/256MiB resource limits, only
after collecting the trajectory replay. No earlier source/default, numerical
diffusion, floor, atmosphere, artificial mass or energy correction is changed.
The full formation sandbox, long-trajectory and local-conservation gates remain
open. `verified_inputs` validates the complete original replay; `constraints`
records native boundary traces, not a replacement for strict solver checks.
