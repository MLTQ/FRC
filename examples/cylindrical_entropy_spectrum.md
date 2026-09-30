# Radial cylindrical acoustic verification

The linearized original PDE is delta_t=−c (r u)_r/r and u_t=−c delta_r,
with u=0 at the outer wall and a regular axis. Its radial frequencies are
c times positive zeros of J1 divided by radius. The dimensionless time is ct/R.

`calculate` assembles the linear Gauss split derivative explicitly; it does not
finite-difference the production RHS to define the matrix. Analytic J1 roots
test accuracy; the positive cylindrical acoustic energy metric tests stability.
Centered differences of the **full nonlinear original driver** at three
amplitudes verify this linearization, preserving their roundoff/truncation limit.

Matrix exponential evolution over two periods is compared with independent RK4
step halvings. All matrices/eigenvalues and raw mode states are retained for
replay. An imaginary linear spectrum is a small-signal component result, not
nonlinear plasma stability or particle confinement.
