# Material magnetic stresses and winding

Source-free, compactly supported full-vector ideal magnetic field carried by
MaterialTwoSpecies. This is a moving magnetic-force component, not a replacement
for finite-mass current dynamics or the existing finite-coil solver. The reference
field must be solenoidal and vanish on the material boundary. Its pullback is
B=F B0/det(F), including the material rotation gradient and cylindrical r/R.
The exterior consequently has exactly zero magnetic field and no matter.
This particular boundary contract excludes penetrating external coil fields.

Magnetic energy is integral V0 |F B0|²/(2 mu0 det(F)). `response` differentiates
this energy analytically at fixed reference flux, returning the meridional force
and angular torque. Both enter the existing canonical fluid equations. There is
no prescribed deformation, fitted force, numerical resistivity or energy repair.
General nodal motion and differential rotation are allowed. `samples` reconstructs
physical B, curl(B)/mu0 and div(B); local computational derivatives use complex
steps through analytic Q1 geometry and the supplied reference-field callable.
The current is a diagnostic of ideal MHD, not an independently evolved finite-mass
carrier state. No unit-test or reporting framework is introduced.

`CompactMagneticField` supplies a reference psi=(Bp/2) r²(1-s²)^4 and
Btheta=Bt*(r/R)*(1-s²)^3 on the actual polygonal reference mesh. Flux gradients
use that mesh's coordinates. Both field and current vanish at the boundary;
J/(en) has a vanishing trace for the n~(1-s²)^(3/2) example. This compatibility
of initial mathematical data does not validate transport/ionization at a cold
plasma edge. Spherical continuum expressions are independently compared by the
trajectory example, including initial polygonal geometry errors.

`MaterialIdealPlasma` adds these forces to the two-species material fluid and
retains its RK4 and strict geometry/species checks. Species entropies are carried
with the bulk in this explicitly ideal single-fluid limit. Full finite-mass
current/Hall/species-relative transport, Joule heat, collisions, exterior field
and finite-coil work remain required for the original goal. Ideal flux freezing
must not be used to claim resistive FRC formation or physical confinement.

Current reconstruction now also solves the physical weak Ampere projection
integral G.J = integral curl(G).B/mu0 with regular nodal vectors
(r*N,0,0), (0,r*N,0), (0,0,N). Boundary terms vanish under this component's
zero-field edge contract. This includes numerical tangential-field face jumps,
which the within-cell curl misses. Both reconstructions are retained and their
difference measured. The projected current is an observable, not an adopted
carrier-current evolution or a guarantee of a vanishing J/n edge trace.
