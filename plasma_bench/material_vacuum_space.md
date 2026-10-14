# Material fields with a penetrating vacuum trace

All three electric, current and magnetic components use compatible material
maps. Poloidal coefficients are nodal psi in the same polynomial x=r0²/2,z0
basis as MaterialMagneticVacuum, with the axis nodal row fixed to zero. Hence
psi is O(r0²) at the axis. Etheta=psi_basis/r0, Br=-psi_z/r0 and Bz=psi_r/r0;
the regular axis limit is evaluated by the polynomial derivative at quadrature.
This replaces the compact r0²*a trace, rather than matching only nodal samples
of two different endpoint interpolants.

`vacuum_trace` selects the actual side/end flux coefficients, in the exterior's
ordering. Normal magnetic flux is continuous over each entire mapped interface.
Toroidal magnetic edge coefficients include boundary values. Physical Btheta=0
in the exterior is a natural weak-Ampere condition, not a strong constraint on
reference Btheta: material rotation can wind a penetrating poloidal field.
Finite-space tangential jumps/current traces remain errors to measure; no
physical singular current sheet is introduced as a separate population.

The original meridional electric/gradient incidence is retained, with zero
toroidal scalar-gradient rows. Actual charge work remains compatible. This space
requires the full vacuum Hessian and Faraday maps in MaterialVacuumCurrent;
using only the interior magnetic Gram matrix would omit exterior flux work.
It supplies no external source coil or formation state on its own.
