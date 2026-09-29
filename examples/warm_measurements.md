# Moving species evolution measurements

Independent Gauss-averaged conservative density/momentum/energy initialize each
case; electron pe starts from exact node samples. Because n is an overlap
average, finite-grid simple-wave U/n is not exactly constant at initialization.
This finite-grid representation error is retained rather than repaired.

Every actual SSPRK2 step contributes diagnostics, stored as compact columns
(including nested limiter columns) without discarding any value. Final arrays retain
both fluid and electron states, exact mesh/volume data and initial states.
Measured nodal pe and cell-average ion errors are explicitly distinguished;
the report independently integrates continuous reconstructed field errors.
Fixed-time controls use exactly their specified count, other cases the fluid
CFL. Any nonphysical stage fails the whole declared study; no hidden floors,
substepping or expected-failure conversion occurs. Source hashes are frozen.
