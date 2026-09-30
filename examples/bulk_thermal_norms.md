# Independent physical moving-thermal measurements

Normal flux primitives are reconstructed as root-product interpolation of the
cumulative face values, independent of production moment solves/integration.
Physical scalar and conservative density bases likewise use independent
root-product formulas. Rate norms include n_t, actual Gamma, both entropy rates,
pressure, temperature rate and local internal-energy rate, with full/first-annulus
L2, absolute RMS and sampled maximum. Weak continuity and physical entropy are
measured separately; neither is a local accuracy oracle.

Trajectory norms use exact characteristics, with the entropy wave rather than
its background as denominator. Particle number, physical log entropy and internal
energy are independently integrated at every saved snapshot. No production
geometry, pressure, energy, density reconstruction or transport routine is used.
