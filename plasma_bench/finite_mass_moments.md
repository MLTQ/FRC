# Exact finite-mass continuum moments

SpeciesMasses retains the exact physical mass fractions and reduced mass;
the explicit electron=0 control is the consistent massless limit. No effective
mass fitting is performed. moments exposes species velocities, drift stress,
bulk/current kinetic energy and full total-energy flux in SI physical vectors.
It includes the mixed bulk/current transport term mu n (v dot w) w and distinct
species enthalpies. It does not evolve fields, close a vacuum or add heat.

reduced_ohm evaluates the exact continuum difference equation from supplied
time/convective derivatives. Callers must include density change in w_t; using
J_t/(ne) alone is wrong when n_t is nonzero. Independent manufactured tests
compare against the two parent species equations and direct species kinetic/
enthalpy flux sums, not a copy of the reduced formulas.
