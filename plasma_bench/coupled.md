# Dynamic gas, prescribed coils and genuine magnetic vacuum

`VacuumCoupledMHD` composes gas-only SSP RK2 half-steps with an implicit midpoint
magnetic/vacuum full-step (Strang splitting). `CoupledState.gas` contains gas
energy only; induced flux spans the larger field grid. `CoupledParameters`
declares gamma, scalar eta, CFL and nonlinear convergence controls.
An optional `max_timestep` explicitly resolves prescribed waveform features;
the fluid CFL cannot detect an arbitrary narrow callback pulse. Waveform/time
refinement remains required. Source filaments inside/on the gas chamber are
rejected because this geometry has no internal solid-conductor exclusion.

External sources are the existing analytic `VacuumField` coil/bias interface.
`step` validates every stage before committing state or work totals. `advance`
retains the accepted state on failure. `diagnostics` accounts for gas plus
induced self-field energy, prescribed-source work, Ohmic heat, Lorentz work,
nonlinear residuals, vacuum-current residual and compatible divergence.

Fluid is confined to a fixed gas-filled cylinder with reflecting slip walls and
ends. Those walls are magnetically transparent: field constraints are solved
through a genuine exterior vacuum. The induced far boundary is finite and fixed.
The exterior's zero-current constraint applies to the induced response; analytic
prescribed source filaments are, of course, external current sources there.
No moving free plasma surface, open gas ends or variable physical transport is
implied. The original experimental Stage 2 gate stays closed until those model
and verification requirements are addressed. This module provides explicit
coupling evidence without relabeling a verification boundary as the experiment.
