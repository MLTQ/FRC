# Recorded continuum diffusion and stiff-time controls

Thirteen stationary-geometry cases compare the actual two-dimensional material/
vacuum induction operator with an independent spherical eigenmode. Records
separate flux, Q1-gradient B, recovered current, local heat, axis and one-sided
boundary current errors, plus full-space energy and heat. Matching, mode,
quadrature, time-cap, tolerance and higher radial-mode controls are independent.
No fluid motion is advanced in this explicitly kinematic component study.

Every diffusion interval uses adaptive midpoint substeps and records accepted/
rejected counts and error/budget estimates. Initial/mid/final field/current/heat
arrays are lossless and hashed. A separate independently diagonalized stiff
low-mode control exposes raw midpoint false reversal despite exact bookkeeping
and measures adaptive tolerance refinement. That control tests time integration,
not a continuum spatial eigenvalue. Source hashes freeze the complete study.
