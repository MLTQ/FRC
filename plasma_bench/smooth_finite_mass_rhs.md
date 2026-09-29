# Smooth coupled finite-mass derivative

Uses smooth_finite_mass_flux for the full coupled bulk mass/momentum/enthalpy
flux and geometric source, then the same physical current induction, drift
stress, current enthalpy and Poynting transport as the frozen baseline. Actual
new mass flux supplies dual n_t and electron entropy transport; no independent
continuity or fixed-density time substep is added.

The short assembly is deliberately separate from the historical frozen RHS so
its original production sources/results remain verifiable. Force/heat validation,
units, regular velocity mapping and FiniteMassRate contract are unchanged.
Electron bulk entropy uses the matching smooth conservative U/n face ratio by
default, with an explicit limited option retained for a controlled comparison.
The lower-order current and electromagnetic representation remain explicit and
must be measured; smooth bulk order alone is not full solver order.
No source callback is present in unforced tests. Inputs never change.
KernelMovingCurrent separates exact curl-free electric backgrounds before its
strict small-signal solve, restores full E for energy flux and checks original
backward errors. The failed original contact solve remains a measured precision
limitation; see kernel_induction for the distinct full-rate replay diagnostic.
SmoothFiniteMassRate extends the previous rate contract with the complete checked
and descriptive induction component diagnostics for scientific records.
