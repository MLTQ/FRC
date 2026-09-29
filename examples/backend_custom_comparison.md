# Existing Hall component comparison

Exercises the unchanged `LinearHallModes` production component with the same
independent cold Hall reference, physical coefficients, helicities and end time
as `warpx_comparison`. Uses both baseline and half time steps. Saved state
coefficients, errors and energy changes demonstrate what our existing component
already covers. This is a four-mode temporal solve without PIC particles or
spatial error; neither runtime nor relative error is a general engine ranking.
It runs locally on the Mac. Production sources are read only and hashed.

Output: `outputs/engine_comparison/custom_hall.json`, used by the comparison
report. This fixture has no plasma-vacuum edge, nonlinear compression or heating.
