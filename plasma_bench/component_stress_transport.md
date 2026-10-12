# Evaluate drift stress by component

`ComponentDriftStress` supplies only the `[..., i, j]` indexing used by the
unchanged `RegularDifferential.stress_rate`. Each requested entry evaluates
`(m_reduced/e**2)*J_i*J_j/n` in that original left-to-right multiplication and
division order. It retains references to the actual sampled current and density;
it does not cache a nine-component tensor or assume floating-point symmetry.
Repeated requests recompute the same scalar field. Other indexing is rejected.

`ComponentStressTransport.pair` copies the scheduled transport method with only
its stress-construction assignment replaced by this adapter. Every force/work
adjoint, actual density, mass, cylindrical connection, face traction, stretching,
advection and self-drift expression is unchanged. There is no quadrature change,
physical approximation, diffusion, filter, floor or repair.

The signed fine replay observed its second transport interval establish a new
peak before entering stress_rate, while the full stress tensor was constructed.
This variant targets that allocation; source-level storage estimates do not
guarantee global peak savings. Original sources remain frozen and defaults
unchanged. Require full coarse/fine old/new physical comparisons and complete
fresh signed replays before refinement, retaining all original strict checks.

The supported contract is intentionally narrow. A future stress_rate accepting
other indexing needs a new adapter; it must not silently fabricate entries.
This is an opt-in instantaneous-work implementation, not a formation model or
stability/confinement result. No unit tests are added or run.
