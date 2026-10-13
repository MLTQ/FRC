# Full spatial contribution diagnostic

The previous goal turn made verified progress: the full failed-trajectory,
tangent and physical-audit replays were recorded and the requested field
display was completed. This continuation keeps the full formation goal open.
All179 preceding Python sources are unchanged at the start. No unit tests.

Use every signed state from the frozen tangent probe: both saved-error and
independent continuum directions at amplitudes1e-4,1e-5,1e-6. Reproduce each
full original RHS and require its fluid, B_t, J_t and E arrays to match the
archived probe exactly before attribution. Record actual input immutability.

Split the native fluid and electric forcing into ten groups:
bulk (including its existing pressure wall reaction), ideal/Lorentz, Hall,
changing-density inertia, current advection, stretching/drift stress,
self drift, electron pressure/entropy, ion pressure/entropy, and the remaining
plasma wall reaction. Save the unconstrained bulk wall value separately.
No component is removed from the solver or used to evolve a modified state.

At each actual signed density the original induction system is linear:
M E = I(n) J_t + f, M J_t = -Fw^T H B_t, B_t = Fw E.
Solve each of the eight non-bulk/non-wall forcing groups with the original
split inverse and1e-10 bound. Restore each species' own compatible potential
to E and check its full Ohm equation. Exact affine species splitting is
recomputed directly; do not subtract large full forces to infer remainders.
Check recombination against the original pairing before using any component.

If a separated solve rejects, retain its input, exception and attempts. Do
not invent its response or weaken bounds. Original full-RHS success and
complete separated-attribution success remain distinct. When all components
return, measure full-array sum differences and the recombined full original
equations. Numerical superposition is a measured property, not assumed.

Physical audit reconstructs each plus/minus contribution independently at
quadrature16/20 and uses their actual physical-rate centered differences.
Keep density, velocity, both specific entropies/pressures, B, J and E in
separate units and norms. Sum signed contributions to check the full
physical RHS. Coordinate/physical projections are signed directional
contributions, not eigenvalues, probabilities or percentages of a cause.
Retain wall momentum reactions in native/physical space even though the
admissible coordinate metric excludes the inactive wall-normal slot.

The known continuum mode supplies independent fluid-term comparisons:
n_t=-div(n0*v), v_t=-grad(dp)/(M*n0)+J cross B0/(M*n0),
s_s,t=-s_s'*(v_z+a_s*J_z/(e*n0)), where a=(-beta,alpha).
dp_s=p_s*(gamma_s*dn/n0+ds_s). These continuum terms are independent of the
production spaces; representation and amplitude errors remain exposed.
The actual saved-error direction requires separate continuum/interface
consistency investigation; this decomposition alone cannot diagnose causality.

Archive every input, component forcing/rate, recombination and full replay.
Use one heavy worker at a time with8GiB peak/256MiB free-disk guards. Preserve
physical masses, boundary spaces, EOS, initial data and original bounds.
The full opposite-sign/time/arithmetic/conservation trajectories remain open.
