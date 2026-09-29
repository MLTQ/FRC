# Current/thermal evolution measurements

Initial and final physical magnetic/current/species states, separate conserved
cell energy and every accepted diagnostic are retained losslessly. The cell
history includes the separate stored-state Ampere/Faraday backward residual,
in addition to the small-signal nonlinear and linear residuals. The cell
energy starts with magnetic, physical current kinetic and both thermal terms,
then changes only by returned face flux. It is never reset from components.
The pressure-unit mismatch, maximum/axis error and inferred conservative ion
pressure are measured separately from evolved positive species pressures.

Every nonlinear failure preserves accepted arrays/time. A nonpositive inferred
ion partition is archived as a rejected trial, not an accepted trajectory.
Failure records contain the actual error and requested/accepted times. No hidden
retry or floor is present. Exact-wave physical B/J and face power errors use
independent direct quadrature; initial/final current coefficients and last-step
electric/face powers permit external physical replay. Completion refers only
to the declared bounded measurement, not FRC formation or a moving plasma.
