# Compatible Hall verification report

Run `python -m examples.verify_compatible_hall` after both measurement protocols
complete. The report first audits all 63 cases and historical baseline identity,
then compares fixed-grid trajectories, nested nonlinear meshes and nonlinear
solver tolerances, and resolves independent wave amplitude/phase projections.
It writes a separate report JSON and three scientific figures.
It never reruns numerical trajectories or rewrites their provenance.

The main engine fingerprint, report source hashes and figure hashes are separate.
`completed` certifies archive/report completion, not physical formation accuracy.
The residuals and denominator conventions retain the first-order physical field
representation and finite mesh differences; small global budgets are not used
as accuracy certificates. Figures must be visually inspected before delivery.
