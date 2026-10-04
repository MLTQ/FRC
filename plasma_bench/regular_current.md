# Current inertia on the actual regular fluid density

`RegularCurrent` specializes the verified variable-current operator's density
representation. It installs `RegularDensity` on the actual regular point state;
it never sends those coefficients through a cell-average reconstruction. All
electric/current sampling, Ampere, inertia action, density forcing and physical
energy routines are inherited unchanged from `VariableCompatibleCurrent`.
Its density gradients therefore return regular point covectors, including the
physical volume weights; density rates have the regular scalar shape.

`electromagnetic_pair` samples the actual reconstructed bulk velocity, Ampere
current and full magnetic field. Its ideal/Hall electric forcing is
−v×B+delta*J×B/(e*n). The physical Lorentz force uses the work adjoint in
`RegularFieldMap`; ideal bulk/magnetic work cancels and Hall work vanishes.
Physical inertia retains the actual supplied n_t. The remaining exchange is
integral mu*J²*n_t/(2*e²*n²). This is explicitly returned, not hidden as heat
or canceled with a fitted force. Complete current advection, drift stress and
species pressure/entropy transport must be coupled before this becomes a
full moving plasma solver. The helper is an instantaneous component, not a
time integrator or a formation model.

Physical masses and density must match the fluid; no adjustable electron mass,
frozen density, temperature reset, diffusion or numerical vacuum is introduced.

`RegularSpecies` now pairs actual species entropy-density rates with additional
pressure forcing on this model. It provides a reduced induction forcing and a
separate compatible electric potential. See REGULAR_SPECIES_VERIFICATION for
its bounded component scope; full moving adoption remains unimplemented.
