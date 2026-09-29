# Compatible magnetic and Ohm matrices

`magnetic_hessian` assembles the fixed compatible energy quadratic on full
periodic psi nodes plus g cells. It uses the same radial mass and exact axial
product integral as compatible_magnetic; boundary psi values remain in the
energy, while the Faraday map fixes their rates. Its sparse product C.T H C
is the positive current stiffness used by the midpoint solve. Tests compare
the matrix action to the independently implemented physical energy gradient.

`transport_matrix` assembles X+R from common physical quadrature. Hall X maps
current coefficients into the weak J cross B/(ne e) electric field; off-diagonal
blocks are installed with exact skew-transposes. R is a positive mass weighted
by nonnegative resistivity. number_density=None explicitly disables Hall, not
carriers. Positive Hall density and finite fields/resistivity are required.
No current floor, anomalous resistance or thermal source is inferred.
