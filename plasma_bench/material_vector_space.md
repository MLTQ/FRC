# Full axisymmetric material field/current coordinates

Purpose: supply the missing poloidal magnetic field and toroidal carrier current
for the moving finite-mass solver. This reuses MaterialCurrentSpace and the
existing compatible polynomials; it is not another verification framework.
Full coupled vector dynamics and complete trajectories are pending.

Reference x=r0²/2,z0 coordinates use the existing electric one-forms Er=r0*Ex*Nz
and Ez=Nx*Ez. The new toroidal one-form is Etheta=r0*Nx*Nz, regular at the axis.
Electric packing is [old meridional coordinates, toroidal nodal coordinates].
The compatible scalar gradient has zero toroidal rows because fields are
axisymmetric. The new curls are Cr=-r0*Nx*dNz and Cz=(2Nx+r0²*dNx)*Nz.

Poloidal magnetic coordinates are a=A_theta/r0, with psi=r0²*a (flux per radian).
Br=-psi_z/r0 and Bz=psi_r/r0. Their joint basis is the curl of the same regular
toroidal one-form, so div B=0 by construction. Toroidal magnetic coordinates
reuse the existing independent two-form basis. Magnetic packing is [Btheta,
poloidal a]. `flux` samples the actual psi associated with those coordinates.

This first full-vector boundary contract is compact: a and its outward normal
derivative vanish at the side and both ends, so all B components vanish there.
The trace null space is computed from actual polynomial value/derivative rows;
row normalization conditions its SVD without changing the physical constraints.
The exterior therefore has exactly zero matter, conductivity and field. This
contract does NOT admit penetrating coil fields or a driven exterior vacuum.
Initialized closed poloidal flux is not a formation result.

`electric`, `magnetic` and `curl` each contain r,theta,z reference quadrature
maps. A full material Piola transform must include rotation and its gradient.
The current solver must retain all magnetic/current stresses and torques, Hall
terms, physical masses and species/heat physics; using these matrices alone
does not establish that coupling or its accuracy.
