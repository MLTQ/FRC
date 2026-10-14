# Moving finite-mass plasma with magnetic vacuum

Extends the full-vector material equations with actual exterior flux, energy,
geometry derivative and traction on the freely moving side and ends. No
exterior fluid/conductivity is introduced. Physical masses, current inertia,
Hall/self-drift, drift stress, both species transport/pressure and Joule/exchange
heat remain in the simultaneous equations. No source coils are included yet.

For boundary selector T, H=H_inside+T.T*Hvac*T and its actual geometry derivative
uses Hvac_dot from the same moving vacuum solve. A toroidal electric coefficient
is the vector potential of the corresponding poloidal magnetic coefficient, so
D[:,Etheta]=-H[:,poloidal] and likewise Ddot=-Hdot. The meridional electric
columns retain the physical interior curl work (exterior Btheta is zero).
This retains both Faraday and the weak Ampere boundary term; adding exterior
energy only to diagnostics would give a different, incorrect current equation.

The same full H,D,Ddot enter finite-current inertia and the unchanged physical
mixed-row residual bounds. The compact parent's interior force is augmented by
minus the exact minimized vacuum-energy coordinate derivative. Its work cancels
the exterior magnetic metric work via `field_response` on the same assembled
exterior geometry. It does not come from the measured residual.
Initial projection includes the target vacuum trace in the magnetic error norm;
evolved states are never fitted to a target or energy total.

The poloidal trace matches everywhere under the actual geometry map. Tangential
matching and current traces are weak finite-space conditions, requiring measured
convergence, particularly as n tends to zero. Toroidal reference boundary values
are free so they can compensate angular winding of penetrating poloidal flux.
Exact vacuum is retained, but heated-edge validity is not thereby established.
The upwind entropy method's measured time sensitivity remains inherited; energy
conservation alone is not a valid local-current accuracy claim.
