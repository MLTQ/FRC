# Conservative MHD state
`MHDState` owns cell conserved variables `(rho,mr,mz,total energy)` and corner
flux. `from_primitive` constructs a state from SI gas variables and the actual
discrete magnetic field. `primitive` recovers pressure without any density or
pressure floors, raising `NonphysicalState` if invalid. `integrated_quantities`
reports exact volume-weighted mass, axial momentum, energy components and face
divergence. Physical/total energy includes magnetic energy exactly once.
The adiabatic index is an explicit solver parameter. No species or transport
closure is inferred by this numerical state module.

The explicit keyword `energy_quadrature` defaults to the historical `cell`
field square and can select positive `staggered` face-field quadrature.
Construction, pressure recovery and diagnostics must receive the same choice,
just as they must share gamma. The state is never converted or repaired after
evolution; pressure recovery only subtracts the declared magnetic energy.
Nonfinite recovered pressure is also invalid. Existing bounded solvers retain
their default cell convention.
The shared quadrature helper owns the magnetic-energy unit conversion.
