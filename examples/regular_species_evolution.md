# Full-transit species evolution measurements

`run_case` advances the actual point entropy through one fastest-axis transit
of the analytic sheared current, retaining physical electron/ion mass fractions.
RK4 uses twelve steps per axial fluid point unless explicitly overridden;
time-refinement cases isolate its error. Endpoint, penultimate state and eight
history samples permit reproducible final-step and full-trajectory replay.
The transit time uses the exact physical constants used in the analytic current.

`audit` independently reconstructs entropy as (S C)/(S n), compares its signal
with the exact sheared translation, and reports axis/wall/max, entropy and
quadratic drift, nonlinear internal-energy drift and initial representation
error. It checks the prescribed current against the analytic curl of Btheta.
Independent quadrature changes are explicit. Internal energy is not an invariant
of the discretized frozen-current subsystem, and there is no energy repair.
