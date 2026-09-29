# Saved-state surface force diagnosis

Restores recorded geometry, velocities, labels, entropy and masses at an explicit
saved time, recomputes the variational force and checks the magnetic solution
against the archived coefficients. Reports equatorial radial acceleration from
pressure, magnetic virtual work and their sum. Separately projects recovered
weak Jtheta times the actual Q1 B to nodal loads with physical volume weights.

These are two discrete magnetic-force constructions, not independent exact
continuum solutions. Their agreement or disagreement helps locate an error;
it cannot prove local momentum convergence, surface closure or formation.
No evolved state, heat, source work or energy is modified. The diagnostic source
belongs in report provenance, separately from frozen dynamic measurements.
Returned profiles and analytic samples use ordinary numeric lists so the
strict JSON record writer can round-trip the complete diagnostics unchanged.

`initial_compatibility` separately derives the initialized FRC's one-sided
surface current and Joule pressure source, checking current by finite Maxwell
derivatives of the independent analytic reference. With continuous B, no other
surface stress/heat terms, and a smooth p=0 material trace, positive edge heating
violates the first-time pressure compatibility condition. This conditional
statement is derived in SURFACE_COMPATIBILITY.md; it is not a general theorem
about weak solutions or a substitute for resolving numerical error.
