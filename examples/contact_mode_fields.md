# Physical reconstruction of the rapid mode

Loads and verifies the exact saved real leading eigenvector from the stationary
360-coordinate contact probe. Differentiates the declared physical fluid
reconstruction: dn=S(dn), dv is the derivative of the actual regular velocity,
and ds=(S(dC)-s_h S(dn))/n_h. The base physical entropy is (S C)/S n. No mode
is fitted, filtered, evolved or substituted into a solver.

The independent root-polynomial sampler reconstructs B and weak Ampere J.
Explicit derivatives of that same magnetic polynomial give its strong curl;
independent current derivatives give strong div J. The exact outer-wall normal
trace is measured without an artificial boundary offset. Reports differences
in normal fluid velocity and toroidal B at the wall, and norm fractions in the
outermost radial field cell, whose boundary aligns exactly with quadrature.
`magnetic_curl` is separately tested on an analytic cylindrical polynomial.
Reports these physical traces
separately from the original passed weak/discrete induction residuals. This
does not omit distributional tangential-B or normal-J interface jumps and claim
a complete Ampere/charge balance; those face contributions still require audit.

The interpolated base pressure is also compared with the analytic stationary
contact. A poorly represented background must not be silently treated as an
exact continuum equilibrium. Quadrature16/20 controls physical integration.
Raw physical perturbations, currents, derivatives and original coefficients are
returned for matched projection and independent inspection. This is diagnostic
evidence about a verification mode, not a laboratory instability prediction.
