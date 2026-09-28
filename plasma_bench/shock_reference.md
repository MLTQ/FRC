# Independent perpendicular magnetic shock

`CylindricalShock` derives a compressive entropy-increasing Rankine–Hugoniot root
for inward flow perpendicular to axial B. Density, pressure, velocity and field
are SI, gamma is explicit. The polynomial follows mass, momentum, total-energy
and induction jumps; no numerical fluid operator constructs the reference.

Away from the standing jump at r_s, each branch obeys r rho v=constant,
B/rho=constant, p/rho^gamma=constant, and H=gamma p/((gamma-1)rho)+v²/2+B²/(mu0 rho)
constant. `branch` solves the scalar Bernoulli relation in log density, selecting
sub-fast downstream or super-fast upstream flow. It rejects radii beyond the
physical choking point. `profile` joins the branches; `flux_invariants` exposes
the independent shock conservation conditions. This is a local steady flowing
solution, not a global cylinder equilibrium and not an FRC configuration.
