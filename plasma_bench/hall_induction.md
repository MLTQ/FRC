# Experimental nonlinear Hall work component

Applies the nonlinear massless-electron cold-electron Hall EMF on the existing
psi/g state through HallCurl. For current j=M^-1(-C.T grad H), sampled J=I j,
raw Eq=(J cross B)/(ne), and e=M^-1 I.T W Eq, the exact assembled identity is
grad(H).T C e=-J.T W Eq=0. This holds for positive spatially varying density,
arbitrary finite amplitudes and all three coupled field components. It includes
both poloidal and toroidal Faraday updates and boundary traces. No Joule heat,
extra fluid force, vacuum density or energy correction is introduced.

The result exposes both rates, coefficients, raw Hall samples, compatible
current, assembled magnetic work and its absolute cancellation scale. Raw
samples and projected electric samples are distinct; only the former are
pointwise perpendicular to J. Conflating them would hide local transfer error.
Density arrays at cells mean a piecewise-constant carrier model; callers can
instead supply explicit common-quadrature values. Zero/negative/nonfinite
carriers fail. psi has the existing periodic conducting-wall/axis convention.

This experimental component is NOT connected to ToroidalMHD or formation.
Independent current/Faraday consistency, local Poynting transport, boundary
behavior and stable time integration are required before that integration.
Exact global Hall work alone is not evidence of spatial or temporal accuracy.
