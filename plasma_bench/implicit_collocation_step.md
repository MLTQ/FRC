# Strict full collocation equation

Solves two real stage increments for the declared solver.a, using actual rounded
physical states and the complete strict coupled RHS. The returned state is the
declared endpoint_stage; solver.weights must equal that row of a (the Radau
constructor supplies this stiffly accurate tableau). The original weighted
endpoint equation is checked independently of the stage iteration.

Uses fixed physical coordinates, unchanged absolute/relative nonlinear and
induction bounds, positive point/sampled states, chord directions with Newton
fallback, residual-decreasing line search and atomic failure. No stage or field
is filtered, clamped or repaired. A complex approximate inverse affects solver
speed only. The equations/tableau determine numerical damping; its physical
accuracy must be measured independently, not inferred from residual convergence.
