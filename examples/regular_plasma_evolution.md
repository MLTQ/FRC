# Complete bounded trajectory comparison

Advances all actual n,q,l,mz,Ce,Ci and magnetic-increment degrees of freedom
together, with physical current inertia and two species pressure/entropy.
There are no manufactured sources. `run_case` covers a full independent wave
period or exact nonlinear contact transit unless an explicit periods override
is recorded. Time steps are fixed and failures propagate; no hidden retries.

`field_norms` reconstructs independently and compares density/pressure
perturbations, actual velocity, B and J. Axis/wall/global/max errors include
initial representation error. Both initial and final norms are saved.
`grid_change` integrates true grid fluid energy and current/magnetic quadrature
energy separately; expm1/log1p evaluate thermal differences without subtracting
large background totals. Particles, species entropies and actual axial/angular
momenta are retained. Wave energy drift is also scaled by the small analytic
perturbation energy, not only a large thermal background.

Raw initial, final, penultimate states, magnetic arrays and dt support replay.
No exact continuous-energy conservation, stability outside the measured cases,
vacuum support or formation physics is inferred. Current degree-five field
reconstruction and RK4 require explicit mesh, quadrature and time controls.
