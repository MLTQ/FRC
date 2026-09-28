# Coupling report and figures

`verify` runs the measurement module, records all configurations and histories,
versions and source fingerprint, then writes JSON, compressed arrays, two PNGs
and `docs/COUPLED_VERIFICATION.md`. The fingerprint covers package Python sources
and both measurement/report modules. Failure aborts; no partial case is counted
as a completed benchmark. Reports from older numerical kernels remain historical
results tied to their own recorded source fingerprint.

`plot` shows independent error refinement, magnetic-domain error, coil work and
energy exchange, and unresolved pulse-flow refinement. A second figure shows
gas density/pressure only within the chamber and current/flux on the separate
magnetic grid; exterior current has no fluid color data. `report` states the
fixed-wall and transport restrictions and the exact response-energy convention.
It explicitly retains the roughly 19% pulse-work change when the finite outer
extent doubles from 2 to 4 m, alongside the separate static Green-function check.
The figures and measured budgets are verification evidence, not FRC formation.
