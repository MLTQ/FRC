# Moving plasma with an evolving vacuum interface flux

Combines the existing bulk motion, physical current inertia, Hall force, current
advection/stress, species pressure/entropy, Coulomb exchange and Joule heating
with VacuumCylinderGeometry and VacuumInduction. The same full weak-Ampere
current enters all terms; the former auxiliary wall-current state is absent.
Density and resistivity are actual evolving fields. Full species-restored Ohm
and induction residuals retain their previous bounds.

The source is a prescribed cylindrical current sheet at C>R. E_source,theta
is -r B_source_dot/2 inside C. The solved electric field is E_ind=E-E_source.
Both E_ind(R) and psi_ind(R) are free. Source work per axial period is
P_source=(2pi L/mu0) B_source mean(psi_ind,t(R))
+ pi C² L B_source B_source_dot/mu0.
It is computed from the source current and linked flux, independently of the
fluid/current/magnetic power sum. Magnetic power includes exterior induced
energy, the nonzero source/induced cross term, and source background energy.
The plasma source work integral J.E_source cancels the B_source_dot part of
the cross derivative. No power discrepancy is repaired or deposited as heat.

Fluid remains impermeable at fixed R and periodic in z. Wall reaction does no
work. The vacuum has exactly no fluid/current/conductivity state. Transport is
still conditional weak-field Z1 Coulomb transport without conduction, viscosity,
neutrals or ionization. This is neither a moving plasma surface nor finite-coil
formation; those parent requirements remain open. End-to-end pulse comparison
must establish physical boundary response before this component is accepted.

Five complete coupled runs and measured spatial/time/model comparisons are in
docs/VACUUM_CYLINDER_RESULTS.md. This bounded periodic coupling is complete;
finite coils/ends, moving free surfaces and formation remain open.
