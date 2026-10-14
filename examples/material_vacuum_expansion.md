# Complete material free-expansion trajectories

Exercises general nodal motion of MaterialTwoSpecies, not an affine-motion
integrator. Initial carrier density nc*(1-s²)^(3/2) and each temperature Tc*(1-s²)
vanish on the moving surface. Entropies, exact particle weights and consistent
physical inertia are retained. No exterior matter or prescribed outer motion.
The initial polygonal geometry differs from the exact sphere/ellipsoid; actual
initial mass/energy and physical-space profile errors are reported explicitly.

Independent continuum reference: gamma=5/3, C=5e*(Te0+Ti0)/(mi+me),
D=a²c/(a0²c0), a''=C*D^(-2/3)/a+Omega0²*a0^4/a³, c''=C*D^(-2/3)/c,
Omega=Omega0*(a0/a)². Sphere with no spin has closed form
R(t)=sqrt(R0²+C*t²). The other axes follow an independent high-accuracy DOP853
solution. These formulas follow the compressible Euler affine physical-vacuum
family described by Sideris (not a plasma transport or confinement theorem):
https://web.math.ucsb.edu/~sideris/pdffiles/explicit-12-09-15.pdf .

Initial radial axis35mm, axial axis35/70mm, nc=1e20/m³, physical Ar+ mass,
Te=Ti=1eV; duration R0/sqrt(C), with optional co-rotation. Analytic total
particles=pi²*a0²*c0*nc/8, U=N*e*(Te0+Ti0)*15/16*D^(-2/3),
K=M*N*(2*a_dot²+c_dot²+2*a²*Omega²)/16. Compare moving boundary, velocity,
density/temperature in material and physical coordinates, both energies,
particle number and axial/angular momenta. The reference positive support part
is never applied to the numerical state. Every accepted canonical state is saved.

This is the zero-current/equal-temperature Euler limit of the intended full
plasma equations, so current/Hall/Joule/exchange terms vanish in the fixture.
It does not implement the missing moving electromagnetic coupling, validate
ionization at a cold dilute edge, or replace the complete finite-coil solver.
No unit tests or new verification framework; actual full free-boundary runs,
frozen sources and4GiB/256MiB guards. Invalid geometry/species reject atomically.

The angle-aware trajectory additionally compares material angular displacement
to phase_dot=Omega in the independent affine ODE. All earlier trajectories and
their pre-angle source closure remain preserved; the fluid equations are unchanged.
