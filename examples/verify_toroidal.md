# Nonlinear toroidal verification report

Checks exact protocol/config/source identities, per-case/collected records,
lossless array hashes, finite saved fields, complete accepted histories and
arrival at the declared end time. Compares time controls only on the same mesh,
normalizing each conservative variable separately by its reference volume RMS.
The finer same-grid run is explicitly a numerical time reference, not an exact
continuum solution. Independent continuum errors remain separate.

Plots exact nonlinear contact/wave convergence, Bessel field/axis/heat errors,
actual time differences, conserved budgets and spurious equilibrium flow. A
second figure compares actual coarse/fine finite-amplitude density and angular
momentum against exact cell averages. Source and image hashes are recorded.
No failed case is silently dropped or relabeled as a completed evolution.
Cell-average norms are explicitly labeled; the diffusion panel also shows
actual reconstructed-field error against physical-volume continuum quadrature.
Transport and wave panels also distinguish cell means from reconstructed
continuum fields with solid/dashed curves and explicit legend labels.
Conservation residual axes use an explicit 10^-15 display scale to avoid
scientific-offset text colliding with the panel title.
