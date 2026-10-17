# Exact coordinate searches within rounded brackets

For finite sorted precise coordinates, float rounding is monotone. Coordinates
whose rounded values differ have the same exact ordering. Two float binary
searches bound the exact insertion index. Vectorized binary search inside each
equal-float bracket compares the original precise coordinates until the exact
left/right index is known. Every tie is resolved, including sub-ULP separations.
No coordinate, field, cell, support or population is changed or discarded.

The owning material model supplies sorted finite edges and finite queries.
Scalar and array query shapes are preserved. Inputs below 512 edges or 128
queries use original NumPy object search to avoid overhead; this cost threshold
does not alter the mathematical result. No tolerance, residual-key construction,
cache or evolving state is present. Sorting/unique remains unchanged.

The prototype retains the entire saved state in the complete 20 ns MPFR argon
run, and Decimal/MPFR physical arrays agree exactly in matched full trajectories.
On one actual full material RK4 step of the accepted 227839-cell state, the
complete output is identical and observed cost is 7.83 versus 12.21 seconds.
These are algorithm-equivalence comparisons under concurrent load, not physical
validation or a controlled performance benchmark. The earlier two-part-key
sorting experiment was slower and remains rejected. Evidence is retained in
`outputs/precise_ordering_cost` and the bracketed argon output directories.
