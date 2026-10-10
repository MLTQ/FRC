# Original physical term and cancellation audit

Reconstructs each archived signed term independently at quadrature16/20 in
bounded radial chunks. Uses the independently tested analytic parent-term
reference. Every term has physical component L2 norms, its own-reference error
when nonzero, error normalized by the corresponding net rate, and separate
axis/wall regions. Zero-reference quadratic terms have no invented relative
own-term error. Electric, raw strong curl and wall-projected Faraday are distinct.

Per-field-cell arrays retain actual, reference and error squared integrals for
all components, plus all pairwise error cross integrals for momentum and
induction. Their sum exposes cancellation instead of adding error magnitudes.
SHA256 streams record every sampled physical term in radial order; chunking
changes memory use only. Full/reduced pressure and its potential are separate;
no curl-kernel information is discarded. The accepted total rate remains an
independent archived comparison. Inertia uses saved numerical Jt, so its error
is feedback from the solution, not an independently isolated forcing cause.

Strong curl excludes distributional interface terms. Previously saved exact
one-sided current/B traces remain necessary. No decomposition changes a state,
operator, residual bound, mass, energy or accepted trajectory. Every raw input,
source closure, radial profile and quadrature metric requires full fresh replay.
