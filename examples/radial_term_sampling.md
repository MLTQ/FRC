# Independent radial chunks for physical term sampling

Slices an existing root-polynomial/Fourier sampler in physical radial rows,
retaining its full axial representation and all fields. This is diagnostic
batching, not a change to a state or spatial approximation. `velocity_rate`
differentiates the actual conservative reconstruction at each signed state.
`negative_curl` uses explicit physical cylindrical derivatives.

`raw_faraday` independently assembles the same difference-first geometric
incidence from coefficient arrays. `projected_faraday` additionally applies the
explicit wall Galerkin radial projection. `wall_projection` constructs that
unchanged production radial mass/trace projection with eight auxiliary axial
cells: its radial matrix is independent of the axial tensor factor. Direct
full-geometry equivalence is tested. Raw and projected curl remain distinct.

Tests compare sliced full sampling, independent physical differentiation and
the original complete discrete operators. No tiny coefficient is filtered and no
interface distribution is silently included in the strong volume curl.
