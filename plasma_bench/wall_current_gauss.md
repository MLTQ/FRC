# Full unsplit time evolution including driven-wall current

Extends the physical state by y(z), the azimuthal wall-current coefficient.
Fluid and magnetic coordinates retain their prior physical scales and exact
boundary-conforming lift; y uses n e cs sqrt((mi+me)/mu), the equal kinetic-energy
current scale. No new residual tolerance is introduced: absolute1e-12 and
relative1e-10 still apply to the full stage equations and endpoint arithmetic.

Both Gauss stages evolve all fluid, magnetic and wall-current variables at their
actual source times. No wall-current equilibrium substitution, splitting,
subcycling, damping, step retry or heat correction is applied. Independent
magnetic coordinates are added before lifting, as in the frozen driven solver.
Continuation requires exact matching of state, B, wall current and time.

The approximate inverse retains the original modal fluid/magnetic blocks and
adds the uniform reference wall-current damping eta*n*e²/mu. It only proposes
nonlinear corrections. Actual changing density, eta, every coupling, full Ohm
and wall-current dynamics remain in each accepted nonlinear residual.
Original secant/Newton fallback and line-search budgets remain unchanged.
Accepted-stage source powers supply integrated work. Full trajectory runs,
not unit tests, evaluate the resulting field/current/heating accuracy.

Measured full24µs160/320-step runs and their remaining spatial/current limits
are recorded in docs/WALL_CURRENT_RESULTS.md. The finer spatial run is still
active; this is not yet an accepted formation or vacuum-boundary model.
