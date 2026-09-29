# Pressure electric and electron enthalpy total-energy flux

`pressure_poynting_flux` integrates the pressure electric field cross B over
primal faces. Ep is poloidal, so only Btheta contributes to poloidal energy
transport. Er=r coefficient and Btheta=r g give the exact annular r² moment;
Ez is at radial nodes. Radial g uses annular-centered interpolation with the
explicit physical wall Btheta, never an implicit copy in place of a trace.

`energy_flux` adds the independently constructed electron enthalpy flux from
`dual_enthalpy`. `energy_rate` takes its conservative primal divergence. It does
not set total energy from magnetic/electron changes. A local difference from
those changes is a numerical ion-energy residual, not physical heating.
The pressure electric gradient is retained in Poynting even when curl free.

Units: face outputs are **integrated watts**, rate W/m³. Axis power is zero;
zero Btheta wall has zero pointwise magnetic/enthalpy power. Nonzero constant
wall Btheta with constant density allows local power but zero integrated input.
This is distinct from the older Hall function's power per unit area convention.
