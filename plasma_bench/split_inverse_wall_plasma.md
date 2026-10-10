# Full original plasma with linear-inverse induction arithmetic

`split_inverse(system)` makes a separate assembly on the same mapping, actual
fluid state/density, geometry, physical masses, boundaries and strict bounds.
Only the induction core electric reconstruction changes arithmetic. Full
species restoration, current transport, physical work and wall reactions stay
inherited unchanged. Cache copies bind the actual current model each call.
No old accepted trajectory or frozen solver source adopts this automatically.
Full captured-input, physical convergence, conservation and atomic controls
remain necessary; this does not establish formation or confinement.
