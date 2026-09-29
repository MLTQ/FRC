# Warm Hall verification protocol

Fifty-eight declared cases: 14 pressure spatial (baseline five grids, three each
gamma1.4/uniform-density/barotropic), six zero-trace diffusion controls, eight
pressure-only evolutions including four time controls, 23 coupled warm waves
(three branches/four grids, opposite axial signs, fixed-grid time, amplitudes,
near-cold electrons), and seven nonlinear coupled mesh/time runs.

Pressure-only evolution holds ions/density fixed to measure the local numerical
ion-energy residual. Coupled cases evolve density, momentum, both temperatures
and all magnetic components. Resistivity is used **only** in the standalone
boundary diffusion reference; coupled warm runs have eta=0. No failures are
declared expected and no case may silently disappear from the collection.
Axisymmetric z-independent diffusion uses four axial cells at every radial
resolution; refinement therefore measures the radial boundary/operator error.
The finest diffusion has a two-half-step time control. Every nonlinear initial
field/pressure/mass coefficient is explicit in the saved `initial` configuration.
Fluid CFL/order/trace-limiter choice, both iteration budgets and the isolated
ion pressure are also explicit; no numerical setting is inferred from a label.
