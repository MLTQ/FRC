# Complete implicit parent trajectories

Initializes the independently analytic contact or finite-mass wave, then advances
the actual full RefinedRegularPlasma state through one full transit/period using
selected implicit midpoint or fourth-order two-stage Gauss. No source, field
correction or parameter fit is used.
The physical coordinate scaling and all solver tolerances are reported.

Saves initial/final/penultimate states, every step's residual/iteration counts,
conservative histories, independent physical initial/final field errors and
strict induction maxima. Grid energy is normalized to the initial complete
fluid/current/magnetic energy and also to wave perturbation energy for the wave fixture. Continuous
energy and quadrature checks remain separate. A progress callback receives
actual accepted history and copied state/magnetic checkpoints; failed trajectories raise and are never called
complete. This initial driver is a candidate until complete controls are archived.
Newton fallback counts and unchanged per-step acceptance limits are also saved;
the maximum residual/limit fraction directly audits nonlinear acceptance.
On failure, the callback receives the last accepted state/time before the rejected
step, including when failure occurs before the first regular progress interval.
