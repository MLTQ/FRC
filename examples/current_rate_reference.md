# Independent continuum toroidal inertia inversion

Uses Fourier z and Chebyshev s=r²/R² collocation to solve for physical gdot:

    L u = u + d0²/N [-4/R² {s u_ss+2u_s-a_s(su_s+u)}
                        -u_zz+a_z u_z],
    N=n/n0, a=log N, d0²=mu/(mu0 n0 e²).

This follows directly from B*=B+curl(mu J/(n e²)), Btheta=r g and physical
Ampere current. Analytic functions of s implement regularity at the axis;
gdot(R)=0 is the physical toroidal trace. No production finite-volume/electric
matrix enters the operator, preconditioner, source or interpolation.

A constant-density Fourier/radial inverse preconditions the variable-density
GMRES solve; its actual collocation residual is checked, not used as continuum
accuracy evidence. Resolution/manufactured-solution studies must establish
that separately. Sources include original Hall, pressure and full nonlinear
inertia. The source's removable axis limit is evaluated at an explicit tiny
analytic probe radius, with an independent probe-size control required; this
is not a simulation density or mesh floor. Cell means integrate the independent
polynomial interpolant in s and Fourier modes exactly in z up to radial Gauss
quadrature, enabling comparison to the physical magnetic rate after inversion.
