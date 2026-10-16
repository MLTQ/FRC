# Controlled common density-partition reduction

`select_partition` removes internal cuts separated by no more than an explicit
physical minimum width. It does not move retained coordinates. Every endpoint
where ANY species occupancy changes is protected, including trace species and
arbitrarily thin supports. Two protected endpoints are never merged. The prior
ULP-only pass remains available and is counted separately.

All species use this one common retained partition. The existing reaction
projection integrates EVERY original Gauss weight and normalized N/P/E moment
into its containing retained cell; no particles, momenta or total translational
energy are discarded. This preserves integrated local reaction charge in the
new space. Density is now a coarser P0 approximation, however, and a pre-existing
nonuniform charge profile and its Gauss field can change. This is explicit
numerical spatial averaging, not physical diffusion or a claim of exact E.

The minimum width controls which nearby cuts may be removed; it is not a bound
on all final widths or on physical error. Protected tiny regions remain tiny.
No positivity theorem follows: U=E-K may still reject after velocity projection.
The width must be refined and local fields/thermal and global budgets measured.

`field_projection_error` compares exact represented old/new P0-charge Gauss
fields, including the empty exterior. New edges are a subset of old edges;
the difference is linear on every old interval, so two-endpoint polynomial
integration gives its exact squared norm. It reports field-energy change,
maximum E change and absolute L2 change. It NEVER changes U or subtracts a
numerical residual from the physical total-energy budget.

Default minimum width0 retains the earlier partition behavior. This optional
reduction alone does not guarantee bounded cell growth, admissible projection,
complete reacting motion, physical model validity or FRC formation.
