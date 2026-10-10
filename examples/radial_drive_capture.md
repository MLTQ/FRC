# Driving terms without accepting a failed induction solve

Evaluates original fluid and electric driving terms at both signed physical
inputs, including cases whose induction solve rejected. No induction solve is
called and no magnetic/current rate is created. All outputs explicitly state
that this protocol accepts no full plasma derivative and preserve the original
full-case status. Inertia is unavailable in a rejected solve and is not filled
with zeros or another case's values.

Verifies exact raw hashes, reconstructed initial fixture and any saved signed
state. When a saved accepted full fluid rate exists, the computed original
momentum/species sum must match it bitwise; absence is explicitly null. Captures
all original bulk/Lorentz/current-stress/species/wall-reaction and ideal/Hall/
density/current-transport/full/reduced pressure coefficient arrays.

These are forcing and momentum consistency diagnostics, not passing total B/J
rates. Same-grid accepted baseline equivalence, independent parent terms,
quadrature and source-frozen full replay are required. The failed original
induction result remains a failure in every table and continuation record.
