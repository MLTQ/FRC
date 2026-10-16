# Compressive numerical bulk viscosity

Optional numerical diffusion for the actual moving material fluid. It is not
physical plasma viscosity, resistivity, ionization or an edge closure. Default
call coefficients are C1=0.5,C2=2; callers must explicitly opt in. No mesh or
thermodynamic state is repaired, and no material is added outside plasma.

With theta=div(v), c=sqrt(sum gamma_s*p_s/rho), and h the smallest physical
quadrilateral altitude, set u=h*max(-theta,0), q=rho*(C1*c*u+C2*u²). Expansion
has exactly q=0. In smooth compression q is O(h)+O(h²); shock diffusion and its
integrated heating require measured convergence. This simple scalar pressure
cannot suppress all shear/hourglass modes or guarantee an untangled mesh.

The nodal force is integral q*div(delta x) dV with cylindrical geometry. Its
power is integral q*theta dV. Positive heat is -q*theta at the same quadrature;
thus mechanical work plus heat cancels without an energy-budget correction.
The caller deposits heat through its actual species representation. The plasma
experiment assigns this numerical bulk heating to ions; it is a declared shock
regularization choice, not a measured physical heat partition. Bulk force acts
on centre-of-mass momentum and does not insert a relative-current drag.

Inputs are one actual MaterialTwoSpecies state and its fields. Outputs include
force, quadrature heating/pressure, power and the unaltered observed work residual.
At a compact edge with rho->0 and finite velocity gradients, q->0. Exterior
matter/conductivity remains exactly zero. Invalid cell altitudes reject.

This uses the artificial-pressure idea of von Neumann and Richtmyer (1950),
https://doi.org/10.1063/1.1699639, with force/heating work pairing of the type
discussed by Caramana, Burton, Shashkov and Whalen (1998),
https://cnls.lanl.gov/~shashkov/papers/comp.pdf. The formula above is the exact
implemented scalar choice; no multidimensional shock accuracy is inferred.
