# Next action: unchanged full RHS under fluid-space refinement

The complete entropy interface budget now reproduces the original operator,
and its fixed-current strong-source error falls from about 55% to 0.029% under
F8→F16→F32 radial fluid refinement at B32/64. This does not prove a cure for the
full coupled growing direction. No exact-mass diagnostic has been adopted.
All preceding physical workers are terminal; do not restart a trajectory.

## Preserve the actual input fields

Use the archived stationary base and both full normalized directions from
trajectory_tangent_probe/context.npz. Keep B32/64, the original magnetic/current
coefficients, 5 mT bias, original physical masses and fixed normalization.
The radially constant base can be replicated bitwise after checking that fact.

Transfer every conservative perturbation component to F16/64 and F32/64 by
exact local polynomial subdivision in x=r²/2, retaining the original Fourier
z coefficients. At parent interfaces, use each child's owning parent trace.
Density, axial momentum and both entropy densities use degree five. Transverse
regular momenta use the original dependent-axis degree-four representation in
the first parent element; omit the new inactive axis coefficient, and retain
the exact normal-wall constraint. Do not insert a nonzero inactive value and
then project the accepted state to hide it. These are diagnostic directions,
not modifications to saved trajectory states.

Independently verify transferred physical density, all velocities, both specific
entropies/pressures, B and current on common q16/20 grids and at one-sided parent
interfaces. Archive coefficient and physical transfer errors. If transfer
cannot preserve the declared reconstruction, diagnose it before evaluating rates.
In particular, constant base density does not justify dropping density-direction
or thermodynamic quotient terms. Keep the original direction normalization;
do not rescale each mesh to a different physical perturbation.

## Complete original operator, with signed controls

For both directions, evaluate the unchanged full split-inverse wall RHS at
both signs of h=1e-4,1e-5,1e-6 on each refined fluid grid, plus the base rate.
Use the original 1e-10 induction checks and all original field/current/species/
wall terms. Save inputs before calls, outputs and any failure immediately.
No relaxed residuals, removed terms, exact-mass substitution or state repair.
Reuse the completed F8 baseline as a comparison; a new full baseline evaluation
must reproduce its raw arrays before it is treated as the same protocol.

Independently reconstruct every physical response, full E and current, and
compare actual differences on the common physical grid, including both species
pressures and conservative-versus-physical differentiation. Keep each component's
norm and signed directional projection separate; these are not certified
eigenvalues. The independent continuum mode remains the separate analytic
control. Retain finite-amplitude and base-stationarity defects. Check original
and continuous power/conservation quantities without global subtraction.

Run one guarded worker at a time (8 GiB/256 MiB), estimate the full raw/replay
payload before starting, collect actual exits, and complete a fresh replay of
all signs/states/metrics. Only then choose complete signed 100 µs trajectories
with same arithmetic, time refinement and full local energy audits. Retain
rejections. A successful scalar response or lower projection is insufficient
to claim coupled accuracy or formation. All original goal gates stay open.
