# Conservative bulk transport for finite-mass plasma

Rusanov transports rho, both linear momenta, total energy and axial angular
momentum. Maxwell stress enters once. The central energy flux contains only
bulk advection of kinetic/current energy and both species enthalpies; full
computed E cross B and drift energy transport are added separately by the RHS.
This avoids substituting ideal E when induction includes inertia and pressure.

Total pressure/velocity/field traces use the existing axis-regular reconstruction.
Separate limited electron pressure fraction and current-energy traces supply
distinct gamma_e/gamma_i and the full physical conservative energy jump.
Scalar limiting contracts reconstructions only; no cell state changes. The
fast MHD speed with gamma p=gamma_e pe+gamma_i pi sets the numerical bulk
dissipation. This is not a bound for every dispersive two-fluid eigenfrequency;
time-step accuracy/stability must also resolve induction and species drift.

Rusanov dissipation and limited reconstruction are explicit numerical diffusion,
not physical viscosity or resistivity. Induction has no added magnetic diffusion.
The resulting species heating/entropy error requires measurement; conservation
alone does not establish a physical shock partition or an entropy inequality.
