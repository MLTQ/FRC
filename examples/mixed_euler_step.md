# Verification-only forcing and independent source-work ledger

Unforced calls use the existing frozen cylindrical RK4 step unchanged.
Manufactured calls add only the analytic continuous forcing at actual RK stage
times. Every stage/final boundary is validated by the production solver;
rejection leaves the input state unchanged.

The physical source work is the independent original EOS/kinetic chain rule
H_U(U_numerical) dot S_exact. RK4 integrates this work and the five conserved
state sources for an **audit ledger only**. Neither ledger enters the fluid
state, temperature or force. Finite-step energy-work remainder remains explicit.
This verification wrapper is not a physical source model for formation runs.
