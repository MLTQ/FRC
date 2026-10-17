# Precise velocity-space extension

Preserve all old independent velocity knots inside the postreaction support.
If a new support endpoint lies outside the old endpoint by no more than its
adjacent old velocity interval, continue that interval's affine basis through
the new edge material. The physical support endpoint remains a fine density
face with exactly integrated N/U and shape-weighted force/work. No old velocity
node or fine density cut is removed. Larger extensions get a new independent
endpoint so extrapolation is explicitly bounded to one adjacent interval.

Fine basis fractions may lie outside [0,1]; partition of unity, exact affine
reproduction, consistent positive mass matrix and matched pressure work remain.
The bounded convex remap cannot use these negative shape weights and must reject
this option. Galerkin and constrained-volume projections can use the extended
basis. Momentum load restriction must account for both shapes at the final
physical endpoint, rather than assigning its entire load to the last node.

This is an opt-in representation of newly occupied material, not coarsening of
existing velocity detail. Inherited internal density cuts/support endpoints are
unchanged. Compare full analytic warm reacting expansion under the current
retained-knot default and both projection methods, then the original charged
20 ns runs with explicit source intervals. Numerical trace/mixing/field/current
accuracy still must be measured.

`inherited_velocity_faces` selects the independent nodes and reports which
support ends are affine extensions. It never removes an old node inside the
new support. Larger extensions retain the normal new-endpoint degree of freedom.
The precise material model must explicitly enable extrapolating shape functions
and use the complete endpoint load restriction; this helper alone changes no
state and is not an adaptive timestep or physical cutoff.

Completed evidence: all six full neutral/charged analytic cases at16/32/64cells
pass in outputs/precise_material_extended, sourcef955cfb7...fe61. Neutral front
errors1.88e-13/1.12e-14/1.42e-15m; finest neutral energy5.68e-14 and population
3.33e-14 relative. Minimum retained width5.87e-30m. A separate16-cell Galerkin
extension control also passes. The charged20ns comparison changes E/J by
1.53e-6/1.46e-4 relative but neutral front by51.34um. Thus the specific analytic
edge-stiffness bug is resolved (bead5.5.7), while general local accuracy remains
open. This helper does not fix clustered interior knots in a newly born species;
that different initialization defect is bead5.5.8.
