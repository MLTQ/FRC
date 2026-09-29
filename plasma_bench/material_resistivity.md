# Physical resistive induction on fixed material geometry

`MaterialResistivity` assembles consistent weighted plasma mass M=integral r³NiNj
and full-space magnetic K on the actual deformed mesh. Vacuum nodes/modes have
zero mass and remain stationary elliptic unknowns; there is no exterior fluid
or resistivity. The explicit scalar eta is in ohm metres, never a transport
closure inferred from temperature. Negative/nonfinite values are rejected.

`induced` completes a material a and uniform bias through the true vacuum.
`current_over_radius` solves M q=(K u)_plasma/mu0 and represents Jtheta=r q;
`current_quadrature` evaluates that current at actual material quadrature points.
The weak residual includes interface-derivative mismatch; current profile and
boundary accuracy need independent refinement, not just small matrix residuals.

`step` uses implicit midpoint on the current frozen geometry with exact bias
increments. It returns a pure `ResistiveUpdate`: new coefficients, midpoint
current, nonnegative cell Ohmic heat, prescribed work, magnetic energies and
independently diagnosed balance/linear/vacuum residuals. Magnetic loss is paired
with 2pi eta dt integral_cell r³ q_mid²; no energy residual is redistributed.
An eta=0 branch completes the changed external field with no heat or label
change; material coefficients are copied directly to avoid cancellation when
subtracting and restoring a much larger bias. Factorizations are reused only
for unchanged geometry/eta and timestep.

The component advances no material coordinates, pressure, clock or accepted
solver state. `ResistiveMovingPlasma` owns composition and atomic commit. Midpoint
is second order and dissipative in energy, but is not L-stable or a monotonicity
guarantee: underresolved high modes can oscillate at large diffusion timesteps.
Accuracy controls and independent current/field tests remain essential.
