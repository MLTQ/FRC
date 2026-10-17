# Compatible material dynamics with explicit precise geometry

Electric force now integrates its exact polynomial on every existing union
interval. For t in [0,1], E=E0+(E1-E0)t and phi_right=a+b*t, b=h/H,
integral(E phi_right)/h=a*(E0+E1)/2+b*(E0+2*E1)/6. The left moment uses
(1-a) and -b. Precise differences precede float coefficients. This is equivalent
to the earlier two-point Gaussian integral up to rounding; no density/field cut
or approximation order changes. Exactly neutral material skips its zero electric
force integral. Full32-cell20ns comparison:144.21s versus233.27s, relative
E difference1.86e-14 and J1.70e-8, with identical986 actual transport steps.
All three full independent charged-affine cases pass; finest relative energy
error3.27e-15. See outputs/reactive_material_argon_constrained_comparison and
precise_material_polynomial (the latter still records the separate unextended
neutral-support failure). Timing is observed runtime, not a hardware benchmark.

`PreciseMaterialSubcells` retains the existing fine P0 density/thermal and nested
P1 velocity equations, consistent mass, physical species scaling, minmod Q,
pressure/work pairing and electric quadrature. It inherits the existing fine
velocity, momentum-load restriction and factored number/energy conversions.
Supports, reference fractions and fine-face
velocities carry50-digit decimal coordinates; width/shape ratios are formed
before conversion to float physical coefficients. No density cut is removed.

`PreciseMaterialState.added` advances geometry within that context. The full RK4
composition in `PreciseMaterialSet.step` uses the same precision for stage sums.
Retaining fine-face velocity differences matters as well as accurate position
addition: rounding those velocities first can erase the represented gradient.
Velocity, thermal state, mass solves, field coefficients and budgets remain
float64. This is not an arbitrary-precision physical-state solution.
Acceleration now solves the same physical mass equation using symmetric
diagonal equilibration from `mass_scaling`, also used by reaction return.
Saved charged-source states exposed acceleration errors as large as1.96e8
relative under the old unbalanced solve. No mass regularization or force
modification is introduced. The existing acoustic generalized eigenvalue
calculation is unchanged; on those states it agrees with its equilibrated
form within7.32e-16 relative.

Optional `velocity_faces` selects an exact nested velocity mesh rather than
rebuilding cosine targets. Its endpoints must span the occupied block and
every knot must already be a precise density face. This preserves a moving
velocity space through source return. Saved `velocity_indices` remain the
authoritative mesh; reconstruct a snapshot with its reference faces at those
indices, rather than assuming the initial target cell count still applies.

`PreciseMaterialSet` keeps the original distinct-support validation and explicit
acoustic/plasma restriction. Union sorting, interval membership and quadrature
use precise positions. Gauss E and its energy use float physical lengths after
exact-coordinate subtraction, including the empty exterior. Sampling and totals
retain the laboratory frame, actual currents and all scaled species inventories.
Negative heat, nonpositive/unrepresentable widths, overlapping same-species
blocks and contact with remote ends still reject. No repaired widths or energy.

This explicit alternative permits continuing precision work while the earlier
scalar-geometry runs keep their frozen source. It is not formation adoption.
The original equations/assumptions are in material_subcells.md and
material_species_set.md. Source projection and saved geometry must use the
corresponding precise path; casting these faces back to float loses the fix.

Explicit extend_velocity permits independent old knots inside the occupied support. Edge cells use the affine continuation of the adjacent old element. All fine cells retain their mass, heat and actual endpoints. Exact shape products, force/work and final-endpoint load restriction include extrapolated shapes; old knots remain independent. Saved state records this option.
# Explicit initial velocity-space spacing

Scalar constants use the selected `precise_geometry.decimal` constructor, so
both supported geometry backends evaluate the same affine basis and forces.

`newborn_spacing_fraction` defaults to zero, preserving the historical nearest
face selection. Nonzero values use `velocity_selection.md` to prevent adjacent
nominal targets selecting almost coincident density faces. This operates only
when there is no supplied velocity mesh; every inherited knot bypasses it.
All density cells remain. The actual initial velocity dimension can decrease
and must be reported. It is an explicit approximation choice awaiting local
convergence, not a density cutoff or a change to physical mass.
