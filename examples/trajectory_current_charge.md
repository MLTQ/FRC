# Physical current divergence and one-sided charge-flux audit

Independently reconstructs both normalized directions from the completed tangent
probe on the actual B32/64 current space. Explicit-cell root-product edge bases
retain the left and right normal traces at every field-cell face, including the
periodic seam, axis and wall. No production trace or current differential is
called. Volume derivatives use the earlier independent TransportSampler.

For each cell, verify the cylindrical divergence theorem: integral div(J) dV
equals outward radial plus axial current flux. These are amperes. Surface
jump integrals have units A²/m² and are separately normalized by R times the
surface integral divided by the volume current-squared integral. Interior
divergence uses R² times its squared integral. The normalization is invariant
under the arbitrary direction amplitude; these are not laboratory currents.

The analytic solenoidal reference Jr=r*k*cos(k*z), Jz=-2*sin(k*z) separately
checks cylindrical face integration against exactly zero cell volume divergence.
It is a quadrature verification fixture, not the physical wall condition.
For actual directions, wall and periodic-seam traces are measured, never set
to zero. Local volume Gauss identities do not prove discrete graph charge,
entropy flux consistency, interface conservation or a numerical instability.

Quadrature16/20, all cell integrals, normal trace arrays, divergence fields and
fingerprints are archived; a fresh full workflow replay is required. Inputs
include the completed/replayed physical contribution audit for cross-checking
its original divergence integrals. No old source/state changes or unit tests.
One worker at a time. Memory/disk are observed at stage boundaries and at the
end (not a continuous guard);8GiB/256MiB bounds still apply. Estimated payload
reserve is128MiB above the minimum free-disk bound before starting.
