# Independent local fallback checks

A manufactured energy conveyor starts with one invalid cell. Removing its
outflow exposes neighboring deficits, requiring four successive replacement
passes. This tests rechecking and propagation rather than only a single
troubled stencil. Integrated mass/energy exchanges cancel, shared CT divergence
remains roundoff, untouched cells retain the exact high-order derivative, and
the original state is unmodified. A weak smooth electric field exercises the
magnetic stencil and diagonal neighbors simultaneously.

A separate concentrated electric curl makes magnetic energy exceed unchanged
fluid total energy. The check must include the proposed corner flux, then
locally replace shared corners without adding thermal energy. A valid proposal
does not evaluate baseline fluxes. Invalid first-order proposals raise without
state repair under both comparison modes. Physical shock/refinement and vacuum
energy evidence are measured separately; these algebraic tests do not establish
accuracy or unconditional positivity of the baseline.

A radial one-dimensional manufactured exchange checks that shared local choices
preserve initially constant Bz/rho. A boundary momentum deficit forces local
replacement of a prescribed outward energy flux; an independent area-weighted
flux sum must equal both the reported outward rate and integrated energy loss.
