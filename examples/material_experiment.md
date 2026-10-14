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
