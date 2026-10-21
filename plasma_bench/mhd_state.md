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

`primitive(..., pressure_override=(mask, pressure))` substitutes declared cell
pressures before the positivity check. Only the dual-energy option of
`ConservativeVacuumMHD` uses it, for cold cells whose pressure comes from the
advected entropy. The conserved state itself is not changed.

## 2.5D state (2026-10-20)
A six-variable conserved state adds L = ρ⟨r²⟩Ω (the cell average of r ρ v_θ)
and g = B_θ/r (the cell average weighted by r dr), with ⟨r²⟩ = (r_lo² + r_hi²)/2.
- `toroidal_energies` returns the cell energies L²/(2ρ⟨r²⟩) and
  ⟨r²⟩g²/(2μ0), which are exact for Ω and g uniform over the cell.
- `primitive` subtracts both from E and returns two more columns, Ω and g.
- `with_toroidal` appends L and g to a four-variable state and adds their
  energy to E.

The four-variable path is unchanged.
