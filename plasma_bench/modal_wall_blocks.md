# Uniform physical preconditioner in axial Fourier blocks

Matching periodic field/fluid point counts make the fixed uniform-background
approximation translation invariant. Exact circulant mass/derivative symbols
come from the original axial matrices. Small radial matrices retain the actual
histopolants, magnetic Hessian, electric masses, physical electron/ion inertia,
wall-admissible lift and fluid radial fluxes. The actual nonlinear RHS is untouched.

For each mode, A=-M^-1 F* H T and K=T* H T+mu/(e²n0) A* M A. Hall uses the
actual radial electric cross mass and coefficient (beta-alpha)B0/(e n0).
Lorentz coupling uses the original fluid/field volume-work adjoint, including
the exact axial trigonometric/polynomial cross integral and real Nyquist mode.
Ideal coupling is the negative kinetic-energy work adjoint in K. Uniform axial
magnetic translation retains the earlier centered preconditioner approximation.

`fluid(mode)` is sparse; `coupled(mode)` returns three small dense coupling/
magnetic blocks. No full-grid dense current or Schur matrix is formed. Every
complex conjugate/FFT convention requires comparison with the independently
derived full uniform operator, including nonzero axial flow and all modes.
This approximation proposes iterative corrections only; it never substitutes
uniform density or averaged species thermodynamics into accepted evolution.
