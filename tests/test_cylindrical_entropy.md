# Cylindrical joint transport verification

Analytic monomial integrals/derivatives verify Gauss quadrature and traces.
Complex differentiation of independently implemented physical cylindrical energy
verifies the actual angular H gradient, boundary inversion and pair jump with
different radii. The r=0 angular boundary is exactly zero.

Original two-dimensional Euler acceleration must converge in the first radial
element; global conserved quantities and **each element's** physical-energy
balance are checked separately. Energy is calculated by the independent EOS
chain rule, with interface energy flux and impermeable boundaries. Density,
axial/angular momentum and both entropy sums conserve; radial momentum is not
a global invariant in cylindrical geometry.

An independent even-polynomial physical reconstruction verifies regular mr/r
and L/r², exact cylindrical volume, and no numerical radius floor. Uniform
axial translation must not acquire a radial or axial acceleration from energy
closure. A rejected nonlinear step preserves its accepted input bitwise.
These tests are component verification, not a formation/stability claim.

Analytic Bessel frequency refinement, physical acoustic energy adjoint, actual
nonlinear RHS differences and fourth-order RK time refinement verify the radial
linearization independently of a mere matrix identity.
