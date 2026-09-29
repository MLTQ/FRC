# Independent toroidal continuum references

`equilibrium`: uniform density, rigid rotation, Btheta proportional to r and
constant axial field, with independently derived quadratic thermal pressure.
Optionally removes the pressure balance for a known J cross B/centrifugal
acceleration. It supplies analytic annular moments, including the axis.

`rotating_contact`: exact finite-amplitude nonlinear solution with rho=f(z-Vt),
Omega=Omega0/sqrt(rho), constant g, Br=Bz=0 and fixed quadratic p(r). Thus
rho Omega² is constant, pressure balances rotation/pinch, entropy and angular
momentum translate at V. Conservative variables are independently quadrature-
averaged, including true continuum energy; no numerical pressure adjustment
hides subcell kinetic-energy error. This is more than a linear wave test.

`torsional_wave`: independent linear Alfven solution, exact through first order
in amplitude with analytically averaged sinusoid and regular vtheta,Btheta~r.
Nonlinear O(amplitude²) feedback is measured separately by amplitude controls.
`toroidal_decay`: weak Btheta=C J1(qr), with J0(qR)=0 giving conducting-wall
Jz=0, positive scalar diffusion and known decay. Fluid feedback is O(C²).
All are verification fixtures, not selected laboratory conditions or formation.

`point_fields` supplies independent continuum density, Omega and g at arbitrary
points, including the regular Bessel axis limit. It allows physical-volume
integration of the actual limited reconstruction error separately from errors
of cell averages. The two norms must not be conflated.
All state constructors accept the explicit adiabatic index; their conservative
energy uses the same gamma as the requested solver configuration.
