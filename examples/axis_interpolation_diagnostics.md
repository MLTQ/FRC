# Physical axis metric and analytic-force diagnosis

On a normalized first element, independently integrate the regular odd velocity
interpolant r*P(r²) through the original Gauss-r point velocities. Its true
kinetic Gram matrix is compared with the diagonal grid-energy matrix using
generalized eigenvalues. A large maximum eigenvalue is an actual physical norm
amplification, not solver time instability. The ratio is invariant under radial
element scaling: refinement cannot remove the matrix mismatch itself.

Apply the full physical work adjoint of the simple analytic force f(r)=r,
then invert that old lumped kinetic metric. The resulting force should be r;
its finite relative error proves that naively attaching the old physical field
quadrature to the point kinetic metric is inconsistent at the axis.

The new positive-node u=vr/r polynomial in x=r² has degree-p−1; its physical
x-weighted metric is integrated exactly by the endpoint quadrature and equals
the positive nodal metric. The same analytic force gives u_t=1. No force fit,
normalization adjustment or energy repair is used. Matrices and polynomial
quadrature are saved independently of either nonlinear force implementation.
This is an obstruction to one proposed bridge, not a claim about all possible
coupling methods or finite-density nonlinear magnetic work.
