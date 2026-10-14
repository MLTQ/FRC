# Complete coupled RHS on transferred fluid directions

Requires the finished and fully replayed `coupled_direction_transfer` protocol.
Re-creates the original split-inverse wall fixture with F16/64 or F32/64 and
unchanged B32/64. The stationary base must reproduce the transferred base bitwise;
both original current coefficient vectors must also match exactly. Original
direction normalization is retained, with no fine-grid renormalization.

Evaluate the base and both signs of h=1e-4,1e-5,1e-6 for both complete archived
directions. Save inputs before the actual full RHS call and check their bytes
afterward, including caught failures. Keep every physical term, original 1e-10
induction bound, physical mass and wall restriction. Return exceptions remain
exceptions and receive no invented rate. A case with a failed minus sign is
recorded without evaluating its paired plus, following the original protocol.

The saved format matches the original trajectory tangent audit: full returned
fluid/B/current rates and E, centered conservative actions, positivity,
constraints and base-stationarity defect. Coordinate projections use the new
grid's declared coordinate metric and therefore are not directly comparable
eigenvalues across grids. The independent physical audit compares components
on the common physical quadrature instead.

No trajectory, accepted-state projection, exact-mass substitution, filtering,
diffusion, numerical atmosphere or changed tolerance. Full fresh protocol replay
and independent q16/20 physical/conservation audit are required. No unit tests.
