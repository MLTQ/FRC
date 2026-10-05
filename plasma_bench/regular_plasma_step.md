# Explicit coupled time-step experiment

`rk4` advances the actual fluid and harmonic magnetic increment together using
four evaluations of RegularPlasma. The bias stays fixed. All intermediate and
returned point states and sampled densities must be positive and finite;
strict induction failures propagate without modifying inputs. No retry, floor,
filter, projection, heat or energy reset is applied.

`field_from_packed` reconstructs the existing bounded harmonic representation;
it does not change boundaries or remove any component. Returned diagnostics
record the worst stage residual and iteration count.

This is classical explicit RK4, not an established time-step selection method
for the stiff moving two-species system. User-supplied steps need actual coupled
spectral/trajectory and refinement evidence. Frozen entropy stability alone
does not justify them. Energy is not claimed exactly conserved in time.
