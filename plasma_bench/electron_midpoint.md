# Electron pressure midpoint

`step` advances the isolated paired pressure/entropy equations by magnetic
implicit midpoint and a temporal logarithmic mean for U. This symmetric
discrete-gradient method is second-order consistent; it is not arithmetic
midpoint for the thermal variable. Picard convergence and the returned-state
residual are checked. Scales include each initial component and its predicted
change, so zero initial magnetic field is not normalized by thermal background.
The physical energy is quadratic in g and linear in U; midpoint preserves its
semidiscrete exchange to the checked solve accuracy. This is not an energy reset.
The finite difference of log U is delta(U)/logarithmic_mean(U_old,U_new).
Using that mean with the entropy-conservative edge flux also preserves the
physical logarithmic electron entropy, to the nonlinear solve accuracy.

Inputs are copied and never mutated. Nonpositive/nonfinite trials or iteration
failure raise. The caller must reduce the time step and establish time/space
accuracy; this simple iteration is not an unconditionally convergent stiff
solver, and centered entropy transport is not positivity preserving. Returned
iteration/residual data belong to accepted trials only. No bulk-fluid state,
temperature floor, transport coefficient or formation boundary is added.
