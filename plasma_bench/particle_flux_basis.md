# Primitive of conservative density reconstruction

`primitive_basis` constructs the same two mirrored degree-3/5, cell-containing
mean polynomials as the existing density map, using independent closed cell
moments in x=r²/2 or periodic z. Integrating each from its left face produces a
continuous normal face-flux lift H. With the conservative face divergence D,
dH=R D and the endpoint fluxes remain the supplied face values. Tensor products
therefore give a flux whose pointwise divergence is the reconstructed cell rate.

This is a kinematic commuting construction, not a pressure or momentum solver.
No matter is placed outside the physical radial domain. Polynomial reconstruction
is unlimited; it has no positivity or shock-stability theorem. Closed moments
and the existing quadrature moment map are mathematically identical but can
differ at floating-point scale; independent tests measure that discrepancy.

`difference` maps normal flux per coordinate area to volume divergence. Periodic
z stores one value per distinct face; radial x has both boundary faces. The
primitive basis accepts degree 3/5 and uniform periodic axial cells.
