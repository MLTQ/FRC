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
The mass solve uses `mass_scaling.solve_mass_banded`: diagonal equilibration
of the same physical matrix and momentum load. This removes the demonstrated
rare-species arithmetic failure without altering density, mass, velocity
resolution or thermal recovery. Negative recovered U still rejects.

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

By default, source return retains every existing velocity knot still inside
the resulting species support, including old endpoints, and adds the new
support endpoints. These knots are protected from optional ULP cleanup. Newly
occupied species/blocks without old knots use the declared initial cosine
resolution. The actual velocity-cell count can grow; no inherited knot is
discarded to meet that initial count. `preserve_velocity_knots=False` retains
the historical rebuild option for explicit comparisons.

This addresses the measured zero-source transfer error: rebuilding cosine
velocity grids changed final current by55–65% in complete20ns controls despite
near conservation. Retaining knots preserves existing velocity-space detail;
it is not a positivity guarantee for newly mixed or newly born material.
Small physical extensions can impose a shorter actual acoustic timestep,
which must be reported rather than bypassed. Local thermal rejection remains.

`last_velocity_knots` records the exact inherited coordinates used by the latest
projection. It and `last_overlap` are cleared before each source stage, so a
failure while gathering/evolving cannot mislabel an older source as rejected.
