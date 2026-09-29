# Linear Hall mode evolution

`LinearHallModes` advances all four magnetic/velocity poloidal/toroidal partners.
It accepts arbitrary finite axial wave numbers and nonnegative radial Laplacian
eigenvalues. This kernel has no Bessel reference or analytic wave frequencies.
Magnetic entries are in T; velocity entries are multiplied by sqrt(mu0 rho).
The curl swaps the two partners with coefficient K. Axial differentiation is
ikz. A positive time mass counts reduced-mass drift energy. The operator's
ideal/Hall part is skew Hermitian and its resistive part negative semidefinite.

`step` uses an increment-form four-by-four implicit midpoint solve for every
mode, leaving caller input untouched and preserving exactly stationary states.
It returns final/midpoint coefficients, physical positive resistive heat,
independent endpoint energy residual and equation backward error. The caller
supplies a validated positive finite physical volume normalization (not an energy correction).
Only the latest time factor is cached. This is second order and A-stable but
not L-stable; large steps can have poor phase/damping despite exact energy.

`hall_cylinder` supplies a numerical radial FE spectrum and Fourier grid.
Parallel reference tests use radial eigenvalue zero explicitly as a Cartesian
limit, not as an axisymmetric cylinder mode. No nonlinear terms, compressible
thermal dynamics or plasma-vacuum boundary are implemented here.
