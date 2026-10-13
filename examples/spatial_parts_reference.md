# Independent fluid-term comparisons

`analytic(mode, r, z, normalization)` evaluates the original linear fluid and
species-entropy terms of the frozen Fourier/Bessel continuum mode. It uses
explicit analytic density/pressure gradients, Bessel derivatives and physical
masses; no production fluid, projection, induction or current operation.
Its summed fluid rates are compared separately to lambda times the mode.
`magnetic_responses` separately solves the original continuum Fourier/Bessel
linear induction system for ideal, Hall and each species-pressure source.
It uses the independently archived curl, physical inertia and analytic pressure
curl; full E restores pressure gradients and each current-inertia response.
The summed B/current/E responses are compared to the frozen independent mode.
This remains the N64 continuum reference, not a new convergence claim.

`interior(sample, state, du, current, bias)` instead differentiates the declared
independent polynomial/Fourier reconstruction of the actual complete error
direction about the saved stationary base. The current is the actual weak
Ampere direction, independently reconstructed with covariant derivatives.
The physical base EOS gradient, including any representation defect, remains.
Full density and entropy quotients are differentiated analytically.

Bulk continuity is -div(n0*dv); momentum includes -grad(dp)/(M*n0) plus
dn*grad(p0)/(M*n0²). Lorentz acceleration is dJ cross B0/(M*n0).
Species entropy is -a_s*dJ dot grad(s0)/(e*n0); bulk entropy uses dv instead.
Pressure rates follow the same physical EOS derivative. All other nonlinear
current transport groups vanish at linear order around v0=J0=0.

For a reconstructed current with nonzero interior divergence, the extra
conservative entropy term is -a_s*s0*div(dJ)/(e*n0); the half-divergence
skew-volume term is retained separately as a diagnostic alternative. Neither
is inserted into the solver or asserted to fix quasi-neutral continuity.

Interior comparisons exclude distributional face/boundary terms and therefore
are not a complete weak consistency or local conservation audit. Their errors
cannot identify a defect without the omitted trace contributions. Current
divergence and base-pressure-gradient norms are explicit, not silently zeroed.
