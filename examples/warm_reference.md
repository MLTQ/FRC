# Independent warm common-velocity continuum solutions

`WarmReference` contains no solver imports. SI argon mass, ne=1e20 m^-3,
Te=5 eV, Ti=1 eV default values define **fully filled ideal verification
fixtures**, not a claim of collisionless/collisional validity in the lab.
Fields are independent of r, have uniform Bz and zero current. Thus these
references have no Hall induction or electron-pressure battery even when
temperatures and density vary along z. The electric pressure gradient needed
by a two-fluid interpretation is curl free; this common-velocity MHD module
does not return that physical electric field.

`contact`: arbitrary translating rho(z-vt), pe(z-vt), constant total pressure
and common velocity; pi=ptotal-pe. Both species obey their exact adiabatic
pressure equation. Finite amplitudes test changing entropy as well as density.

`simple_wave`: one Riemann invariant is fixed. v=a sin(k xi),
c=c0+(gamma-1)v/2, rho/rho0=(c/c0)^(2/(gamma-1)), ps/ps0=(rho/rho0)^gamma,
with z=xi+[c0+(gamma+1)v(xi)/2]t. Newton inversion is checked and rejects
characteristic crossing, |t a k (gamma+1)/2|>=1. This is an exact finite-amplitude
compressing Euler solution before shock formation, not a linearized wave.
The constant axial field changes neither dynamics nor reference equations.

`cell_means` uses independent Gauss quadrature, including conservative momentum
and total energy; it never constructs energy from average primitive fields.
Separate 12/18-order comparison estimates quadrature error. Reference PDE tests
check continuity, momentum, total energy and both species pressure equations.
