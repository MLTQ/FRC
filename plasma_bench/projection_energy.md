# Energy paired with momentum projection

For the same mass matrix `M`, source load `b` and projected velocity `v`, write
`K_c = 1/2 integral_c rho_target |v|²` and
`W_c = integral_c p_source · v`. The proposed thermal projection is

`U_c = E_source,c - W_c + K_c`.

Consequently `delta E_c = 2 K_c - W_c` and
`sum(delta E_c) = vᵀ(Mv-b)`. The exact Galerkin solve therefore conserves global
energy. This redistributes local total energy consistently with momentum;
it does not preserve each original fine cell's total energy separately.

When source density is constant within a retained fine cell, the same formula
is `integral(U_source + rho |u_source-v|²/2)`: source thermal energy plus
unresolved velocity variance. It is nonnegative and invariant under a uniform
velocity shift in that limit. General source-density variation, floating-point
cancellation and a changed density representation can spoil that guarantee;
negative output must still reject. No temperature floor or global energy reset.

The added variance is **numerical mixing**, not physical atomic heating.
The caller must record its thermal gain and local energy redistribution, assess
local errors and convergence, and retain the original cell-total-energy option
for explicit comparisons. Finishing a run or preserving global energy is not
sufficient for adopting this research alternative in formation calculations.

`projected_frame_energy` consumes fine-cell energies, endpoint momentum loads,
target endpoint velocities and actual target kinetic energies. It returns the
new thermal array plus signed local total-energy changes in the same units.
It is applied to every cell as a defined projection, not only to failed cells.

Background on the related staggered-remap conservation/positivity problem:
[Kenamond, Kuzmin and Shashkov (2021)](https://doi.org/10.1016/j.jcp.2021.110254).
The formula above follows from this code's Galerkin orthogonality; it is not
an implementation of that paper's intersection-distribution algorithm.

Handoff evidence (2026-10-16): the explicit candidate succeeds on the saved
retained-knot10ns rejected source, but Ar9+ reaches2.882045e6eV from~0.199eV.
Global energy error3.03e-15 relative therefore does not justify adoption.
See outputs/reactive_material_argon_retained_knots_10ns/projected_frame_local_errors.json
and docs/GOAL_STATUS.md. No full candidate trajectory has been run;
local_total remains the default. Beads5.5.4/.5 track the local velocity and
energy defects before full E2E refinement under5.5.6.
