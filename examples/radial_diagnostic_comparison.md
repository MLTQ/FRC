# Driving baseline, quadrature and rejection controls

Requires every common saved coefficient between the accepted full-term capture
and the independent driving capture to be bitwise equal; lists exact common
keys. Compares all component actual/reference/error L2 norms across independent
quadrature16/20 and records the worst absolute and relative cases with their
scales, including tiny zero-near terms. No universal relative claim is made.

Requires the instrumented rejected input to reproduce the original first inner
failure exactly. Full failed-case status is retained. Freezes comparison source
and every input hash; a complete rerun through wall_diagnostic_replay compares
all non-runtime metrics. No state, term, residual or acceptance threshold changes.
