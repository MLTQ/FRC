# Actual full coupled stationary spectrum

Constructs the complete centered finite-difference Jacobian of RegularPlasma at
uniform physical density, unequal species temperatures, rest and uniform axial
bias. All active fluid and magnetic increment degrees of freedom are included;
the inactive axis slots and fixed flux nodes are excluded. No branch is filtered
or fitted. Epsilon controls are required because this is a numerical Jacobian.

The metric is independently differentiated physical total energy: the local
fluid EOS Hessian, physical regular-coordinate kinetic Hessian, plus magnetic
H and Ampere^T I Ampere. Scaling uses n, sound speed and equal-energy magnetic
units. A Cholesky transform measures the symmetric part and all eigenvalues
in this positive metric. Raw matrices and spectra are saved.

The reported 2 sqrt(2)/spectral_radius applies to RK4 on an imaginary spectrum.
It is a diagnostic for this exact frozen mesh/state, not an automatic bound for
nonlinear moving density, a finer mesh or formation. Finite-difference real parts
and energy defects must be retained and controlled, not silently projected out.
