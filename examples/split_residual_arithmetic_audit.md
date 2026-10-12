# Independent arithmetic audit of the complete captured rejection

Audits every recorded iterate from a source/input-verified and matched original
rejection. Uses Decimal.from_float at precision 2200 to compute the exact sum
of the finite stored binary64 electric terms minus their rounded sum, then
compares it with the recorded two-sum low part. This is a full actual-field
audit, not a synthetic unit test; all coefficients participate.

Reports signed projections and vector remainders for the captured residual
decomposition, keeping norms distinct from attribution. It validates saved
vector identities and records all file/source hashes. Complete postprocessing
replay is required. It neither accepts a rejected rate nor changes any solver,
threshold, physical equation or state. An explained arithmetic gap is not a
verified numerical repair or a claim of physical formation accuracy.
