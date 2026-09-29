# Continuous contact norms

Independently integrates the actual continuous representations of electron and
ion pressure and temperature over physical cylindrical volume. Electron pressure
and density are Q1 in x=r²/2,z; ion pressure and density are cellwise constant.
Temperatures are evaluated as pressure/(e number density) in eV. These norms
include interpolation/projection error and differ from the nodal electron and
cell-mean ion metrics in `finite_mass_contact_case`.

`measure` consumes a saved finite-mass state and its evaluated component split.
Physical interpolation and volume quadrature are independent of production maps.
Returns relative L2, dimensional RMS and sampled maximum for pressures and
temperatures, plus absolute B−B0 and J drift norms since their targets are zero.
The audit uses both initial/final states and quadrature orders 10/14. No state is
changed and there is no relative division by a zero magnetic/current target.
