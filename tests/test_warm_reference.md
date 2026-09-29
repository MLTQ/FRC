# Independent warm continuum reference checks

Fourth-order finite derivatives check the independently evaluated contact and
finite-amplitude simple wave against continuity, Euler momentum, total energy
and **both** species pressure equations at two gammas. Exact energy quadrature
is compared at orders 12 and 18; characteristic crossing is rejected. Actual
solver evolution must refine toward these electron-pressure solutions rather
than merely preserve an implementation-specific invariant.

The initialization helper is for tests only. Measurement initialization remains
in its own provenance-tracked module and uses the same independent reference
fields, not this test's thresholds as a substitute for scientific records.
