# Conservative MHD state
`MHDState` owns cell conserved variables `(rho,mr,mz,total energy)` and corner
flux. `from_primitive` constructs a state from SI gas variables and the actual
discrete magnetic field. `primitive` recovers pressure without any density or
pressure floors, raising `NonphysicalState` if invalid. `integrated_quantities`
reports exact volume-weighted mass, axial momentum, energy components and face
divergence. Physical/total energy includes magnetic energy exactly once.
The adiabatic index is an explicit solver parameter. No species or transport
closure is inferred by this numerical state module.
