# Slow-wave continuity diagnosis

Separates the actual limited Rusanov mass flux into central and dissipative
parts, against independently averaged linear parent-species continuity. The
diffusive norm is its magnitude relative to the exact physical rate; the other
norms measure error. Their scalar norms are not additive.

Compares a separately tested cubic annular momentum interpolant using both the
actual finite-amplitude initial state and exact linear momentum averages. This
distinguishes spatial error from nonlinear amplitude truncation. Saves complete
rate arrays at two densities and five meshes; it never advances or changes a
physical state. The coupled solver does not use the cubic diagnostic.

Checks the independent manufactured heat's regular-axis limit at several small
positive radii, densities and phases. This is a reference evaluation control,
not a production density/radius floor. Full enthalpy, momentum and kinetic work
must be consistently reconstructed before a higher-order mass flux is adopted.
