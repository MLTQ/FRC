# Alternative evaluation of the same raw and projected curl

`evaluate` applies the existing assembled wall curl, sparse projection after
raw curl, dense radial projection after raw curl, and a subtraction-first
raw stencil. All operations consume full actual captured electric arrays.
No variant changes a state, solver, residual or acceptance rule.

`exact_sparse_action` uses exact binary64-to-Decimal conversion and 2200-digit
products/sums, then rounds each row to binary64. This suffices for finite
binary64 product exponents and the short raw-curl rows. It is an independent
reference for the stored raw sparse coefficients, not an exact projected
operator or continuum derivative. Dense projection remains binary64.

The subtraction-first stencil uses the same rounded reciprocal mesh factors
as the original CSR matrix; its summation ordering differs. Projection follows
the unchanged magnetic-mass-orthogonal wall map. Full captured-state audits and
fresh full replays verify the workflow; no unit tests or solver adoption claim.
