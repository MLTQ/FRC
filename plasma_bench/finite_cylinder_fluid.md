# Finite-cylinder fluid with impermeable end caps

FiniteCylinderFluid retains the regular state (n,n vr/r,n vtheta/r,n vz,Ce,Ci),
physical mass, two entropy equations, radial SBP elements, dependent axis trace
and impermeable radial wall. The axial coordinate is now a finite Gauss-Lobatto
interval including both ends. It is a global polynomial, not a wrapped Fourier
or periodic finite-difference field. The existing symmetric two-point physical
flux is differentiated with Q=W D and the explicit end boundary terms.

At each cap the numerical flux is pressure only: no particles, advected momentum
or entropy leave. The two endpoint axial-momentum values are inactive and exactly
zero. Their reactions do no work. The end pressure force and additional discrete
constraint reaction are reported separately; the latter must decrease for a
resolved smooth solution. Total axial momentum is not conserved in a chamber
whose end pressures differ. Particle number and zero wall energy flux retain
the physical SBP balance. There is no added viscosity, filtering or temperature
reset; invalid density/state values reject.

The axial pressure derivative uses pressure differences to preserve a large
uniform background at floating precision. The inherited radial equations are
unchanged; their periodic term is explicitly disabled. The inherited RK4 step
advances all species and momenta together and preserves the inactive slots.
The intended first chamber is sealed with impermeable end caps. This is not an
open-end loss model, sheath model or moving free plasma boundary.

This is the fluid portion of frc-vacuum.2. Magnetic/current maps and the open
vacuum must also acquire finite axial traces before the full coupled solver can
use it. The periodic field map cannot be substituted on these new axial nodes.
Independent full acoustic evolution measures end pressure response, convergence
and conservation before the finite-coil coupling is connected.
