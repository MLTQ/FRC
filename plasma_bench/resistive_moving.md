# Resistive moving material plasma

`ResistiveMovingPlasma` composes physical scalar-resistivity induction/Ohmic heat
with `MovingPlasma` material dynamics. Its explicit eta>=0 has SI units ohm m.
It is a prescribed coefficient, not a species/temperature transport model.
The eta=0 path uses the original ideal step without changing labels or entropy.

For eta>0 each accepted step is resistive half-step, full ideal Verlet step,
resistive half-step. Each resistive stage uses declared error-controlled midpoint
substeps (default local tolerance 1e-4); counts, rejected trials and error estimates
are reported. Extended time advances only in the ideal stage; resistive
stages hold the bias at that stage's time. This symmetric nonautonomous Strang
composition avoids advancing the prescribed source twice. Source work remains
in the ideal stage. Cell heat changes kappa=p V^gamma by
(gamma-1) Q_cell V^(gamma-1) at fixed geometry; no heat is assigned to vacuum.

All stages run on copied mutable arrays and scalar accounting. Only complete
success commits state, flux labels, entropy, time, steps, source work and heat.
Failure after either half-step leaves the accepted state intact. Fixed cell and
nodal masses and immutable geometry mappings are shared. Diagnostics retain
physical cumulative cell heat and separate resistive substep residuals from
the full moving/source response-energy error, which is not forced exact.

The ideal component retains its declared interpolation/axis and matching errors;
raw implicit midpoint diffusion can oscillate for underresolved fast modes;
step doubling controls field/current/heat error but is not a positivity theorem. This
solver has no shock capture, remapping, arbitrary spatial coil source, Hall,
ionization or thermal-transport closure. Independent diffusion/current/local
heat, temporal-splitting and dynamic source/refinement checks are required.
