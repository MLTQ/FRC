# Full plasma assembly with balanced induction refinement

`balanced(system)` makes a separate assembly on the same physical mapping,
fluid, geometry, masses, boundaries, tolerances and cache policy. Only the
induction refinement implementation differs. All full species/current/work
coupling, potential restoration and original Ohm check are inherited unchanged.
No state or accepted-rate projection, physical diffusion, filter or repair.

Operator caches contain only geometry; shallow copies bind actual density at
each evaluation. Original UniformWallPlasma and old archives remain unchanged.
Full captured-input and independent known-solution/equation controls are
required before using this opt-in assembly in further rate experiments.
