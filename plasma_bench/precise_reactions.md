# Moving reaction transfer on precise supports

`PreciseReactionOverlap` keeps the inherited actual atomic-moment evolution and
physical observations while gathering on decimal-coordinate union intervals.
Points and interval membership are resolved before float physical coefficients
are formed. Empty exterior rows remain exactly zero. Physical masses, charges,
donor moments, initially empty products and source normalization retain their
existing checks; no population is discarded to simplify the geometry.

`project_precise_reactions` integrates all original quadrature weights into
common piecewise-constant species amounts and energies, solves the consistent
P1 momentum projection, and recovers each actual U=E-K. Negative thermal states
reject. All source fields remain float64 with per-species scaled N/U; precise
positions and nested fractions carry the distinct support geometry.

Only the explicitly requested earlier ULP cleanup is available. Its tolerance
is coalesce_ulps times the float64 spacing at the remote-domain scale, declared
separately from50-digit geometry. Every species occupancy endpoint is protected.
There is no physical minimum-width coarsening option. Cleanup can still alter
represented local fields; old/new Gauss and field energy differences are observed
without compensating heat or subtracting the residual from physical energy.

`PreciseReactiveMaterial` inherits the original source-half / material-full /
source-half and complete budget contracts. It replaces gather/return only.
It preserves actual variable populations and supports throughout motion rather
than keeping a frozen reaction fixture. General charged accuracy, mesh growth,
source dynamic range and all original formation gates remain separate issues.
