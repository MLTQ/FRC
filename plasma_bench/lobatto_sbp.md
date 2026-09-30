# Endpoint-collocated SBP and dependent axis trace

Degree-p Gauss-Lobatto quadrature has p+1 points, both endpoints, and exact
polynomial integration through degree2p−1. The barycentric derivative gives
Q=W D and Q+Qᵀ=diag(−1,0,...,+1). p must be at least two.
Derivative roots are checked real before discarding zero imaginary dtype parts.

The positive x nodes supply a separate degree-p−1 regular radial/angular
velocity polynomial. Its trace at x=0 is dependent on those positive values.
There is no independent axis velocity with zero kinetic quadrature weight.
This distinction is necessary for a positive physical velocity norm and must
be retained by any driver using `axis_trace`.
