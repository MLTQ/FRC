# Independent mode on separate meshes

Generalizes the frozen continuum-wall fixture to an explicit spec with separate
fluid radial elements/axial points and field radial/axial cells. Physical fields,
stationary contact, physical masses and independent N=64 continuum mode are
unchanged. The fixture samples conservative differential fluid variables and
analytic poloidal flux/cell-averaged toroidal field, with the declared one-time
wall-conforming initial approximation. No numerical eigenmode is fitted.

Returns `(mode, fluid, base_state, system, bias_field, du, db)`. `du/db` are unit
scaled *linear* perturbations. Both signs with an explicitly reported small
epsilon are physical candidate states; nonlinear trajectories must keep the
actual conservative variables and initial approximation. Mesh compatibility,
fixed traces and positivity validation remain in the original constructors.
Existing frozen fixture code and raw sources are unchanged.
