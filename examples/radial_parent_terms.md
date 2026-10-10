# Independent physical term reference

Evaluates original linearized pressure and Lorentz accelerations, plus ideal,
Hall, species-pressure and physical current-inertia electric fields and their
negative cylindrical curls. Uses only analytic background profiles, original
Fourier coefficients and Bessel identities. Imports no production momentum,
pressure projection, current, induction or time solver.

For F=(p(z)J1,t(z)J1,w(z)J0), -curl F=(t'J1,-(p'+q w)J1,-q t J0).
The pressure term explicitly retains the perturbed density multiplying the
background pressure gradient. Its negative curl is q sum a_s(p_s' dn-n' dp_s)
J1/(e n²), not the curl of a fitted equilibrium. Inertia is mu lambda J/(e² n).
All terms have physical units. Their large cancelling magnitudes must be shown
alongside the small net rate when interpreting a normalized error.

`terms(mode,r,z)` returns dictionaries `electric`, `induction`, `momentum`.
Independent finite-difference curls, full-field reconstruction and original
mode equations verify signs, geometry, density variation and cancellation.
It retains the source reference's Nyquist convention; use the resolved N128
archive for radial16, not the unresolved N64 control.
