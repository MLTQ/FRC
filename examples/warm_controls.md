# Supplemental warm controls

Six declared runs refine the strong contact and Mach .12 simple wave to
nz128/256 and compare nr4/8 at exactly the same 512 time steps. The latter
separates radial dependence from the base radial control's changed fluid CFL.
Stronger-gradient refinement tests the physical electron entropy defect exposed
by the main protocol; it does not relabel that defect as exact entropy work.

Uses the unchanged frozen main engine/reference/measurement sources. This
supplemental protocol/driver has its own SHA256 in each record and collection;
both identities are checked during assembly. No base records are overwritten.
