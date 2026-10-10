# Linear mass-inverse electric reconstruction

Opt-in balanced core preserving the original operator, residual definitions,
strict limits, physical masses, boundary space and all iteration budgets.
Computes E=M^-1 f + M^-1(I Jt) instead of M^-1(f+I Jt). These are the same
linear equation, but separately rounded sparse solves need not distribute
accurately enough when a small curl is reconstructed from a larger E.

The forcing inverse is cached once for this solve and used identically in the
original RHS and physical E reconstruction. Every original Ohm/Ampere/operator/
dual and physical Faraday check still accepts or rejects normally. No residual
is omitted or replaced, and no bound is relaxed. This is an arithmetic proposal,
not a continuum accuracy fix, state projection or modified physical forcing.
Old balanced-core sources and rejected records remain unchanged. Independent
analytic tests, exact captured-state probe and full replay are mandatory.
