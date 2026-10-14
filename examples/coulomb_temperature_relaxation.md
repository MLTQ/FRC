# Complete electron-ion relaxation run

Evolve a uniform, stationary, fully ionized Z=1 argon fixture in the actual
moving finite-mass equations. Inputs: n=1e20 m^-3, Te=5 eV, Ti=1 eV,
mi=39.948 atomic mass units, both gamma=5/3, R=0.035 m, periodic length0.24 m,
F4/16 and B8/16. Bias, initial current and prescribed resistivity are zero.
Only the new physical temperature-exchange source is nonzero analytically;
the full fluid, current, pressure and magnetic equations still run. No
numerical damping or state correction is introduced.

The exact uniform solution has n constant, v=J=B=0 and Te+Ti=6 eV. Independently
separate the scalar temperature ODE and integrate
t = integral(Te(t),5) [nu_E(T,6-T)*(2T-6)]^-1 dT. Adaptive scalar quadrature
and a root solve in log(Te-3) supply each reference temperature without a
near-equilibrium quadrature singularity; no production collision
routine or full-system time integration supplies that baseline. The temperatures
stay below10 eV, so the low-Te logarithm branch applies throughout. Compare the
initial production rate separately with the NRL rounded cgs/eV coefficient;
its rounding difference is not a discretization error.

Run400 us with40/80 Gauss steps to measure time convergence of the actual
thermal evolution and energy drift. Save every accepted state/stage and solver
history. Independently reconstruct fields, temperatures and physical energy at
q20; report spontaneous motion/current/field and entropy change. This uniform
fixture has no spatial structure to establish spatial convergence, drift or
magnetized transport validity. It is not an experimental comparison or formation.

Original acceptance tolerances apply. Discrete own-RSS/free-disk checks after
setup and each accepted step impose1 GiB/256 MiB; there is no continuous guard
during native calls. Sources/runtime and results use the existing record helper.
No unit tests, new verification framework or automatic replay is added.
