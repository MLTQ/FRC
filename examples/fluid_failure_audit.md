# Last-accepted-state and actual rejection audit

Requires the frozen calculation identity and raw hashes. Reconstructs every
failed case's actual retained last accepted state, replays its stored history,
and independently integrates its physical energy at orders18/22. Repeats the actual
failed RK4 step and checks the original error plus unchanged accepted arrays.
The audit runs on the original platform; its own source is archived separately
without changing the frozen calculation closure.
The initial audit rejected a remote failed state's physical energy comparison:
the unstable state's quadrature error is now retained explicitly as a warning,
including higher-order density positivity and quadrature convergence. This does
not relax a solver or hide the discrepancy behind a larger passing tolerance.
Rejection/history replay and physical quadrature accuracy have distinct results.
