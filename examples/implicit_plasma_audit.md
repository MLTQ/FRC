# Independent physical implicit-trajectory audit

Reconstructs the complete initial/final density, velocity, separate species
pressures/temperatures, J, B and full E at independent quadrature orders 16/20.
Continuum contact/wave norms retain axis/wall/max diagnostics. Integrates the
original physical bulk, both internal, current and magnetic energy, including
the constant-bias cross term. Thermal differences use stable conservative-change
reconstruction and expm1/log1p, without changing any accepted state or energy.
The audit treats neutral laboratory fill settings as unrelated to these analytic
fully ionized verification fixtures.

Reports continuous physical energy drift separately from grid energy and,
for waves, perturbation energy. Initial total energy includes harmonic bias.
Particle and both species-entropy changes are normalized by initial particle
count; point/sampled density and both temperature extrema remain explicit.
Integrated energy does not establish the unresolved local face budget.

Requires original source/runtime and raw-array hashes, freezes its own complete
closure, and refuses overwrite. Optional full replay requires every saved array
and every non-runtime trajectory metric to match bitwise. A replay or audit
exception leaves no terminal-success summary. No replay is implied when the
full_replay field is false. The separately archived chord-only failures remain
unchanged and are not mislabeled successful trajectories.
