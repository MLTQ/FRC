# Nonsingular-coordinate physical reconstruction

Density, axial momentum and both entropy densities are degree-p polynomials
in x=r²/2 through actual Lobatto point values. u=vr/r and Omega=vtheta/r are
polynomials through active nodal velocities q/n and l/n: degree-p−1 in the
first element with its dependent axis trace, degree-p elsewhere. Reconstructed
physical momenta are n*r*u and L=n*r²*Omega. The axial velocity is mz/n.
All point values at active DOFs are retained, not fitted, clipped or filtered.

Independent root-product bases and periodic Fourier resampling evaluate this
declared representation. It is regular at exactly r=0; positive/finite density
is checked at sample points. Pressure and higher quadrature must be checked by
the caller. Continuous energy/angular integrals can still differ from grid
quadrature and must be measured, not declared exactly conserved.
