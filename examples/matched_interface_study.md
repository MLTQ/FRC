# Frozen physical interface-consistency study

Evaluates the actual complete bounded-plasma RHS at independently prescribed
continuum perturbations on separately declared fluid/field meshes. Each case
has volume and exact one-sided physical face budgets at quadrature16/20,
including current divergence/jumps, particles, both entropy densities and
original total energy. Every term, cell rate/balance, face integral and original
RHS coefficient is saved, with expanded physical-face fingerprints.

Source/runtime/companions freeze before the calculation; existing outputs are
never overwritten and terminal success requires the entire explicit case list.
Replay must compare every coefficient/integral/fingerprint and non-runtime
metric bitwise. Large physical cell residuals are reported, never corrected by
subtracting global work or assuming pointwise charge continuity. These are
instantaneous consistency measurements, not a changed solver, full trajectory
local proof, moving vacuum boundary or formation experiment.
