# Full original physical cell and interface budgets

Uses independent volume and exact one-sided face quadrature on aligned field/
fluid meshes. Integrates actual full-RHS particle, species entropy and total
physical-energy rates, then original parent-species fluxes through each cell's
own interior traces. Current continuity has no evolved charge variable;
its measured cell outflux, strong divergence and normal jumps stay explicit.

Reports all cell residuals, full strong-current divergence, outer-wall flux,
periodic/internal left-minus-right jumps and pointwise absolute jump integrals.
The periodic seam is internal; axis area is exactly zero. No physical trace is
averaged or forced continuous. Geometric summation identities and the independent
current divergence theorem test the bookkeeping separately from physics.

For particle/entropy/energy, relative cell L2 uses the combined volume-rate and
outflux-density norm. For current it uses the independently reconstructed full
physical current-gradient Frobenius norm; a zero charge rate never normalizes
itself to a spurious constant-one result. Absolute units, maxima, integrated
powers and all cell arrays remain. Expanded face fluxes have exact fingerprints.
These measurements diagnose local consistency; they do not repair fields,
energy, interface fluxes or the underlying weak equations.
