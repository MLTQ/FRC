# hybrid3d_coils.py

Builds a compatible discrete magnetostatic basis for the known external
coils. Analytic/cubic coil sampling followed by a discrete curl can leave
spurious current inside the gas, particularly next to a nearby winding.

`VacuumCoilBasis(box, gas)` solves B_source=−grad χ inside gas, with the
original normal flux on its boundary. Only source fields are constrained;
the evolving plasma response is never averaged, filtered or projected here.
For the present poloidal sources, annular gas has no axial Bθ circulation.

The correction δB is divergence-free and confined to gas. A periodic discrete
Poisson inverse gives δA=(−Δ)⁻¹curl δB; then curl(A+δA)=B_source to roundoff.
Thus the source vector potential, imposed induction, plasma-current
subtraction and circuit reciprocity remain one compatible field basis.
The artificial periodic boundary lies inside the outer conductor padding.

`hybrid3d_sources` and `hybrid3d_circuit` use this when gas geometry is supplied.
This source solve changes no ion, temperature or plasma mode and does not
assume that the actual plasma is current-free. It enforces Maxwell's vacuum
condition on the independently known external field only.
