# Directional physical-error and arithmetic report

Creates a six-panel standalone PNG/PDF from terminal frozen measurements.
Separates complete accepted rate errors from driving-only diagnostics on the
preserved failed input. Full term norms, errors, radial squared-error fractions
and independently reconstructed floating-point residual gaps retain their
different normalizations. Inertia includes numerical-current response feedback.
The stored radial coordinates are field-cell centers; radius normalization
recovers the outer face from the uniform center grid.

Requires accepted full rates for the four plotted mesh controls; refuses to
plot a rejected derivative as accepted. The doubled radial-field case uses the
explicit opt-in split-inverse arithmetic with all original checks. No claim
of full trajectory accuracy, exact vacuum, formation or confinement is made.
Uses the driving audit's explicit full-reference-rate normalization. Freezes
source, hashes all plotted inputs and output files, refuses overwrite.
Visually inspect the generated figure; the terminal record alone is not a
layout check. Failed plotting leaves no successful terminal summary.
