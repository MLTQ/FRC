# Exact separable electric-mass inverse candidate

Keeps every original CompatibleGeometry basis, quadrature point/weight, mass
matrix, Hessian, incidence curl, boundary trace and sampler. Only mass-inverse
evaluation changes: if M=Gr kron Gz and f reshapes to F, solve
Gr X Gz^T=F with one radial and one axial LU. All modes and coefficients remain.
No field fit, truncation, physical-mass change or approximate preconditioner.

TensorMassInverse.solve requires a complete finite real vector in the original
row-major electric ordering. The original assembled mass is retained for all
physical equations and strict residuals. Rounded assembled Kronecker entries
and separated triangular solves have different finite-precision paths; bitwise
mass-inverse equivalence is not assumed. Acceptance limits remain unchanged.

TensorMassGeometry reproduces the original initializer expressions and inherits
all sampling/packing/gradient methods. StrongWallGeometry can wrap it normally.
The old geometry is immutable and remains the default. The new class is opt-in
until full original-versus-new physical workflows, matrix/input identities,
independent fields/work and convergence, and complete fresh replays establish
adoption. A coarse complete pilot is not fine-grid accuracy certification.

This addresses costly global two-dimensional mass factorizations, not all
memory use: assembled sparse matrices, current metric and nonlinear quadrature
arrays still exist. Record actual resources; do not promise a peak-memory bound.
