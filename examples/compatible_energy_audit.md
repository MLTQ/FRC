# Independent finite-step energy decomposition

For the constant-bias shock, response energy is a quadratic function of fluid
total flux plus a linear function of conserved fluid energy. The Hessian uses
full-space vacuum completion minus fluid magnetic energy. For stage rates k_i,
the exact RK defect is dt² times the Gram contraction with
`b_i b_j - b_i a_ij - b_j a_ji`, plus the weighted semidiscrete defect.

`quadratic_step_defect` independently completes each rate through magnetic
vacuum, forms these bilinear products and predicts the energy error. It does
not read the measured energy change or alter the solver. Time-dependent sources
are explicitly rejected: their affine completion requires additional terms.

`run` observes the original accepted RK stages without changing their fluxes,
times, limiter reports or budgets. It checks the prediction against each actual
response increment and persists cumulative results at the standard four shock
outputs. Incomplete records are marked as such. This diagnoses time integration
versus a spatial mismatch; it does not repair energy or establish physical
accuracy. Source/package/audit fingerprints and final diagnostics are saved.
