# Identical initial-state comparison

`from_split_initial` converts an unadvanced split fixture's gas energy to
total cell MHD energy by adding its existing staggered-field cell energy.
Geometry, sources, induced vacuum field, gamma, physical eta, CFL and source
cap remain identical. Advanced states and optional added shock diffusion are
rejected: this is initial-fixture construction, not an energy repair or a way
to inherit a split solver's evolved budgets. Reuses independent shock, smooth
advection/diffusion and source fixtures for controlled solver comparisons.
Conservative-only interface controls use their declared defaults, since the
split fixture has no corresponding numerical interface-flux parameter.
Explicit parameter overrides support first-order, limiter-disabled, interface
and temporal control runs; they are validated at construction and recorded.
