# Exact retained-knot projection diagnosis

Focused study of the actual saved rejected source in
outputs/reactive_material_argon_retained_knots_10ns. Reuses all original decimal
velocity knots; no cosine-mesh reconstruction, coarsening or altered source.

Separates (1) the production balanced solve, (2) 100-digit tridiagonal elimination
of the rounded production matrix/load, and (3) independent 100-digit assembly
from raw quadrature moments, precise physical coordinates and target P0 density.
The independent assembly uses analytic P1 product integrals per fine cell and
direct source-point shape integration. It does not call production load
restriction or mass assembly. Decimal arithmetic uses 100 significant digits.
Neutral 0 and Ar9+ are compared explicitly, including a constant-velocity
mass-row-sum arithmetic probe. That probe alone is not physical verification.

Three snapshot transfers preserve full states: Galerkin/projected-frame,
bounded paired moments, and energy-constrained volume velocity. Reports include
all species' local velocity/temperature ranges, source-point velocity error,
normalized number/momentum preservation, global energy and numerical mixing.
The two candidate transfers are research comparisons, not adopted physics.
Complete time evolution and independent limits are evaluated with the existing
reactive_material_subcycling, reaction_remap_control and precise_material tools.

Run python -m examples.reaction_projection_diagnosis FRESH_OUTPUT. Saves source
provenance, raw states, measured JSON and scientific PNG/PDF. No unit tests,
provider input redistribution or generic verification framework is introduced.
