# Dynamic plasma/vacuum regressions

Tests enforce actual exterior magnetic-only degrees of freedom, positive heat
partition, uniform-state balance and closed-end gas conservation. Independent
Sod, polynomial J cross B, translating compact flux and broadening Gaussian
solutions verify hydrodynamics, magnetic force, advection and diffusion.

The magnetic subsystem is checked using independently evaluated quadratic field
energy, gas/internal energy, prescribed-source work and exterior current. Both
ideal and resistive cases are covered. A complete two-coil pulse must accelerate
and heat gas while conserving the measured response budget. Nonlinear failure
must leave state, time and cumulative work untouched.
The explicit source-resolution time cap is enforced, and source filaments in
the gas-filled chamber are rejected because no internal conductor is modeled.
Both fixed and infinite multipole exteriors are covered by the ideal/resistive
magnetic work and internal-energy checks, including energy beyond the grid.

An infinite-space circular-loop integral separates finite-box truncation from
spatial error using exact boundary data in a verification-only Poisson problem.
No test treats energy conservation, zero discrete divergence or a passing smoke
case as proof of formation, free-interface validity or magnetic-shock accuracy.
