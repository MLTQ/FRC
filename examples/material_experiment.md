# Run a configured moving-plasma experiment

Usage: `python -m examples.material_experiment CONFIG OUTPUT`. The JSON describes
the tube, pre-ionized inventory, initial compact plasma, temperatures, grid,
prescribed resistance, bias and pulse. See `plasma_bench/material_experiment.md`
for equations and assumptions. The supplied metre-scale argon case is exploratory,
not a recommended gas, fill, plasma preparation or hardware specification.

The complete coupled system advances all actual states with inherited RK4. Every
accepted state, resolved inputs, exact input config and frozen source are saved.
Both source-work integrals and Joule heat are independent trapezoidal observations;
they never change the simulated state. Source self-energy baseline is excluded
as documented in MaterialCoilCurrent; mutual I*lambda remains in its response
budget. Failures retain all accepted states and identify the attempted step.

Outputs include actual field/current components, bulk expansion, temperatures,
flux extrema and the plasma volume fraction where Bz opposes the **local source**
Bz. That sign statistic is not an O-point/X-point search, closed-flux certificate,
axis measurement, FRC formation or confinement evidence. The pulse itself may
reverse the source field. All sampled minima/maxima use plasma quadrature points.

Existing conditional regime diagnostics are evaluated at the actual maximum
electron-drift/thermal-speed sample. Their length/time are explicitly supplied
reference scales, not automatically resolved gradient scales or validity proof.
Other important samples and nonlocal/neutral/kinetic physics can remain invalid.
Small Coulomb logs, the cold edge and prescribed resistance stay visible.

This runner adds experiment inputs and actual simulations, no unit-test suite,
replay framework or tolerance-based acceptance gate. It stops at wall contact
instead of silently imposing a new boundary condition. Source modules must stay
unchanged while its live trajectory is running.

Independent source groups are recorded with their actual currents, centre-field
contributions, induced linkage and both separately integrated work contributions.
Totals use the full summed physical source field. Legacy scalar current/linkage
values are null for multiple independent drives; no fictitious common current is
invented. The midplane boundary radius and particle-weighted RMS radius are also
observed, distinguishing a local pinch, bulk extent and the maximum surface radius.
None of these radii proves particle confinement or models escaping particles.

With optional integrator="split", source-response work and its group integrals
are null because instantaneous upstream pressure at zero current does not define
that split work integral. Direct external work -sum(Idot*lambda), actual Joule
heat and group external work remain independently integrated. Existing split
diagnostics are recorded; all accepted states and any failure remain preserved.

The optional declared numerical bulk viscosity is observed separately: cumulative
numerical ion heat, peak numerical pressure and actual mechanical/thermal power
residual. Its heat uses the same independent trapezoidal observation as Joule heat;
neither measured integral repairs the state or closes the budget by construction.
The config and resolved metadata retain the coefficients; omitted means no bulk
viscosity. See the complete wall-shock baseline in `material_wall_shock.md`.

For paired collisions the historical `joule_heat_j` means collisional relative
flow converted into total species heat, not prescribed eta*J². Independently
integrated `collision_heat_species_j` retains the actual signed electron/ion
heats, including thermal exchange, while numerical ion heat stays separate.
Records also retain collision entropy production, force/heat power mismatch and
the maximum local relative-relaxation rate times dt. These observations neither
repair energy nor certify a full-system timestep limit. The whole coupled source
work balance remains the trajectory conservation measurement.
