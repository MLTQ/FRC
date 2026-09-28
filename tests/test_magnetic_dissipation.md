# Numerical diffusion and midpoint-equivalence checks

Uses analytic regular cylindrical compression to require the exact divergence,
O(h²) numerical diffusivity, and zero resistance in expansion, rest and vacuum.
Compares direct and cached-vacuum implicit solves with physical/numerical
resistivity independently enabled, under both exterior options. Independently
integrated field/gas/internal-energy changes check source work and the separate
physical/numerical heat partition. No test accepts a repaired pressure or energy
state, and no equivalence claim extends to an unconverged midpoint iterate.
