# Conservative density at compatible field quadrature

`CompatibleDensity` reconstructs **actual number-density cell averages** with
the existing degree-3/5 annular and periodic moment polynomials. Each candidate
stencil contains the evaluation cell, so the reconstruction preserves its mean
when quadrature resolves the polynomial. This is not a nodal density, a second
continuity evolution, or a replacement for the caller's mass flux. The caller
supplies n_t=rho_t/(mi+me) in exactly the same cell-average representation.

`sample` applies the fixed linear reconstruction R to any signed rate.
`adjoint` returns R.T W f with physical cylindrical volume weights, without
division by cell volume. Thus sum(adjoint(f)*n_t)=integral f R(n_t) dV. This
pullback is required when assigning density-dependent current-energy work.

Cell means and quadrature samples are copied and immutable. Nonpositive means
or reconstructed samples reject; no floors or fallback are used. Positivity is
checked at quadrature points, not proven between them. The smooth polynomial
map is unlimited and is not a shock or vacuum reconstruction. Geometry and its
quadrature/degree are fixed for the object's lifetime.
Finite input coefficients whose reconstructed values overflow also reject;
the same check applies to signed density rates.
