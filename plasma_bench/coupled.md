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
through a genuine exterior vacuum. The default comparison fixes induced flux
at a distant boundary; the multipole option below instead extends to infinity.
The exterior's zero-current constraint applies to the induced response; analytic
prescribed source filaments are, of course, external current sources there.
No moving free plasma surface, open gas ends or variable physical transport is
implied. The original experimental Stage 2 gate stays closed until those model
and verification requirements are addressed. This module provides explicit
coupling evidence without relabeling a verification boundary as the experiment.

`vacuum_boundary="multipole"` represents the induced field to infinity through
decaying modes coupled by magnetic energy. `multipole_modes` controls truncation;
moving the matching surface and refining modes/mesh are separate accuracy checks.
The fixed zero-flux option remains the default comparison for historical fixtures.
`beyond_grid_induced_energy_j` accounts for exterior energy outside that surface;
the modal stationarity residual is separate from physical nodal vacuum current.
An infinite magnetic exterior does not remove the fixed gas chamber restriction.
`midpoint_method` selects the original direct linearized solve or a Picard update
with cached vacuum completion. Both solve the same implicit equations and obey
the same nonlinear residual/atomic-commit contract; comparisons measure their
field, gas and budget agreement as well as runtime. Picard is the default, with
a tighter default relative tolerance of 1e-13; direct remains an explicit reference.
`shock_diffusion` is an explicit, dimensionless numerical coefficient, zero in
the original centered comparison. Numerical magnetic heating and maximum added
eta are reported separately from physical resistivity/Ohmic heat. Enabling it
does not constitute a transport closure or validate physical reconnection.
The time-step accuracy cap includes the maximum current numerical resistivity,
as well as physical resistivity, keeping its fixed-point iteration controlled.
