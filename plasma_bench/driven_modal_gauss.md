# Full Gauss stages with prescribed source time

`DrivenModalGauss` extends the frozen LiftedModalGauss step to a nonautonomous
axial-flux model. It requires `step(..., when=t)` and checks the input field
against Bbar(t). Every Gauss stage evaluates its actual source time t+c_k dt,
with that time's complete field and inductive forcing. Even identical stage
state guesses must evaluate both times. The endpoint uses Bbar(t+dt).

Independent magnetic coordinates, physical variable scales, full nonlinear
and induction bounds, positivity checks, original line search and all fluid
arithmetic remain. The modal inverse's initial background is only an
approximate preconditioner; it is not the time-dependent physical field.
No splitting, fixed-temperature replacement, state repair or step retry.

Each residual returns the source powers belonging to its own stage rates.
Accepted power values are returned as diagnostic stage_source_powers_w and
source_work_j=dt*(P1+P2)/2. Rejected trial/Jacobian evaluations cannot replace
these accepted powers. The original Newton helper only consumes the first
tuple member and is unchanged. Physical energy/work error remains measurable.

Continuation also requires an exactly matching source time. Carry the returned
endpoint time forward with the same arithmetic as the field to avoid comparing
an independently rounded source time. Coordinate memory commits only after
the original endpoint checks pass. Frozen autonomous integrators remain intact.

Use with the explicit AxialFluxDrive models. Filled-cylinder flux boundaries
and current-trace limitations are in axial_flux_drive.md; this does not provide
finite external coils, exterior vacuum or physical ends. Verification is an
actual driven evolution against the independent diffusion reference.

Three full driven runs are recorded in docs/DRIVEN_FLUX_RESULTS.md. Fixed-eta
Bz improves on8→16 cells, but space-time current error remains28.06%→21.14%
and heating7.87%→4.47% low. The evolving-eta case is a physical control
comparison. The wall-current approximation remains an unresolved solver issue.
