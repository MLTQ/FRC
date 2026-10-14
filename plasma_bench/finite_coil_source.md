# Finite circular-coil source

An explicit list (radius, axial center, signed ampere-turn multiplier) shares one
prescribed waveform in amperes. Analytic elliptic loop fields and Atheta=psi/r
are evaluated at actual common physical quadrature. The affine magnetic field
keeps induced coefficients separate from the exact source field; source current
is never mislabeled as uniform axial magnetic bias.

linkage(J)=integral_plasma J.Aunit dV is the reciprocal induced flux linked by
the source per common ampere. The response energy includes I*linkage and the
source response power is I*d(linkage)/dt. Only the divergent prescribed filament
self-energy baseline is excluded. Source E=-Idot*Aunit drives physical induction.
The volume reciprocity expression and the vacuum field evaluated at coil
locations converge to the same linkage; finite-space discrepancies are physical
representation errors to measure, never compensating heat or repaired work.

There is no conducting or massive exterior. The analytic source filaments have
prescribed current; plasma-induced vacuum field remains the solved source-free
exterior. Finite wire geometry and hardware circuit dynamics are not modeled.

Seven full24us coil pulses now exercise this coupled path. Independent limits,
spatial/time comparisons, transport flags and remaining exterior-resolution
limits are recorded in docs/FINITE_CYLINDER_RESULTS.md.
