# Continuation acceptance and atomic rejection

For all three continuation variants (including bounded proposal work), two actual full plasma steps are independently
substituted into the original Gauss stage and endpoint equations. The second
step exercises the previous-step proposal; unchanged input arrays are required.
A deliberately exhausted changed-dt step must leave both physical input state
and every previously committed continuation-memory value unchanged. These are
solver-equation and transactional tests, not growing-mode convergence evidence.
