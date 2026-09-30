# Independent cylindrical physical reference and integration

`initial` samples the independent analytic reference into actual point n, n vr,
n r vtheta, n vz and species entropy densities. `physical` and `rates` use the
original EOS and cylindrical kinetic energy, independently of production field
or energy-gradient routines. Complex input permits independent differentiation.

`regular_sample` declares a finite physical reconstruction between nodes:
piecewise polynomials in r² for n, mz, C_e, C_i, mr/r and L/r², and periodic
trigonometric interpolation in z. Thus mr is odd and L vanishes quadratically
at the exact axis; no radius floor or divergent L/r field is allowed. Root
products construct the polynomial basis independently of the SBP derivative.
The reconstruction passes through all evolved point values. It is a diagnostic
physical representation, **not the polynomial in r used by SBP differentiation**.
Regular interpolation does not prove that evolved states approximate a smooth
axis solution. Near-axis original-equation errors remain separately required.

`measure` evaluates independent high-order physical energy and field errors.
Nodal Gauss conservation does not imply conservation of this continuous energy.
Both their difference and successive higher-quadrature differences must be saved.
Positive nodal values do not guarantee positive interpolated density/pressure;
invalid sampled physical states fail rather than being floored or clipped.
Velocity error is RMS divided by the physical background sound speed, including
when the analytic velocity is identically zero. This avoids undefined relative
errors at standing-wave turning points. Other fields use their actual L2 norm.
