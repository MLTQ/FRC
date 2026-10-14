# Compatible toroidal-field/poloidal-current material space

Uses the existing finite nodal/histopolant polynomials in x=r0²/2 and z0.
Electric one-forms have Er=r0 Ex Nz and Ez=Nx Ez. Toroidal magnetic two-forms
have Btheta=r0 Ex Ez. The raw Faraday incidence R maps electric coefficients
to magnetic rates and annihilates compatible nodal gradients. All end traces
are finite; no periodic wrap is present.

Btheta=0 is imposed at the side and both ends by the exact reference polynomial
trace null space Z. This persists under a meridional Piola map. Axis regularity
is supplied by r0. There are no poloidal B coefficients in this invariant sector.
A regular, source-free, axisymmetric toroidal magnetic vacuum has Btheta=0;
no exterior fluid or conductivity is introduced. This sector excludes coil
fields and cannot describe FRC formation by itself.

The geometry holds sparse quadrature maps for E, B and actual curl(E).
Moving physical Galerkin projection is assembled by the current operator;
the old fixed-reference boundary projection is not reused after deformation.
Reference matrices only are cached. Current/induction uses actual deformation,
particle weights and density at every stage.
