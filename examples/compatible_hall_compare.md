# Archived trajectory comparisons

`state_difference` compares compatible magnetic fields, density, total energy
and the three conservative momentum components. Its denominator is the fine
reference's **change** from its initial state, avoiding dilution by the uniform
magnetic bias, density or pressure. A zero reference change returns null.
Angular momentum is converted with the cell's exact mean r² before its norm.

`restrict` integrates all nested fine cells with annular volume weights and
samples coincident corner flux. This is a conservative/corner comparison, not
an independently exact continuum field or pressure projection. No evolved state
is altered. `comparisons` collects fixed-grid time controls and nonlinear mesh
comparisons against the saved nr64 run, with actual step counts.

`wave_modes` constructs independent exact g cell-average templates at zero and
quarter-period phase, then projects the evolved toroidal coefficient in the
physical r²-weighted magnetic norm. It reports fundamental amplitude, wrapped
phase error and residual outside that two-template mode. These separate damping
and dispersion from the larger magnetic reconstruction error. Zero-Hall cases
with identically absent toroidal perturbation are explicitly excluded.

The report auditor owns archive identity and source validation. This module is
report-only, separately fingerprinted; it does not relabel numerical provenance.
