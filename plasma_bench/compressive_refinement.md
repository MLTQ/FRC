# Refinement of unresolved material compression

`refine_compressive_cells` is an explicit research option for a
MaterialSpeciesSet of full-cell MaterialSpecies blocks. It does not accept the
scaled or separate-velocity-space subclasses used by the reacting solver.
It adds cells/nodes without merging, removing populations or changing supports.

Use the existing minmod residual compression jump:
C=max(0,max(-dv,0)-h*min(max(-g_left,0),max(-g,0),max(-g_right,0))), g=dv/h.
Select C>threshold*c, with the physical sound speed from each actual thermal
state. Affine compression has zero residual; expansion is not selected. The
threshold is a numerical resolution parameter, not viscosity or physical heat.
Refinement is allowed only if both children meet the explicit minimum child
width. Existing narrower cells remain. Width-limited requests are reported.
Exceeding the explicit maximum cell count fails instead of silently truncating
the selection. No mass change, sound-speed relaxation, temperature floor or retry.

For selected cell[a,b], add its midpoint, halve N and U, and linearly interpolate
nodal three-velocity. The represented P0 density/thermal and P1 velocity fields
are unchanged in exact arithmetic. Particle number, thermal energy, momentum
and the exact integral of kinetic energy are conserved. Gauss field and empty
exterior are unchanged too. Before/after total budgets measure rounding, without
compensating it. Cell widths/amounts that cannot represent the split reject.

Refinement changes subsequent spatial dynamics and numerical viscosity. In
particular, preserving a coarse P0 field does not recover unresolved gradients
or guarantee accuracy. This option needs full independent shock comparisons
and equal-time plasma-current/field refinement measurements. It is not enabled
in formation or reacting transport. No unit suite is its validation criterion.

The first full20ns constant-child plasma experiment preserves energy but creates
current oscillations:128 initial cells grow to3664 electron cells without
improving the approximately34% difference from the1024-cell reference. Retain
that failure as evidence; more cells alone do not imply improved transport.

`thermal_reconstruction='limited_linear'` is a separate explicit candidate.
Compute parent pressure p=(gamma-1)U/h and minmod the neighboring secants
between physical cell centers. Boundary slopes are zero without ghost matter.
For slope s, assign Uleft/right=U/2-/+s*h^2/[8(gamma-1)]. This gives child
pressures p-/+s*h/4; minmod bounds them within neighboring pressure extrema
and preserves positivity for nonnegative input on unequal cells. Actual invalid
thermal states still reject without clipping. Parent total U is conserved,
but the local P0 thermal field deliberately changes. Report the absolute
thermal redistribution and maximum child-temperature change at every event;
this is numerical reconstruction, not physical heat. N, P1 velocity, support,
kinetic energy and Gauss field retain the constant-child split contract.

Completed full plasma evaluation rejects this variant too: 2,904 final electron
cells retain a 32.77% current difference from the 1,024-cell reference, with
large new oscillations. The constant case gives 34.17%. Both remain research
options only. See docs/COMPRESSIVE_REFINEMENT.md for independent shock results,
full-duration budgets, redistribution and the new fixed 2,048-cell comparison.
