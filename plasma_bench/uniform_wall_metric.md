# Uniform physical inertia metric as an inner preconditioner

The actual induction operator is K(n)=H+A* I(n) A. The old magnetic-Hessian
preconditioner omits current inertia; this module approximates the inverse with
K(nref), using a fixed physical volume-mean density and the same physical masses.
Radial matrices and exact axial circulant symbols produce small SPD Fourier
blocks in the exact wall trace-null space. All modes are retained.

The original normalized Krylov coordinates are handled by positive diagonal
similarity factors. Small modal roundoff is Hermitian symmetrized only inside
the approximate inverse; original accepted equations and strict residuals are
untouched. For n in [a*nref,b*nref], positivity bounds K(n) between
min(1,1/b)K(nref) and max(1,1/a)K(nref), motivating this preconditioner without
freezing the actual evolving density. No density floor or mass adjustment.

`solve` maps a real normalized Krylov dual to an approximate correction. Its
SPD symmetry/positive work and uniform-state exactness require independent full
physical-operator tests. The reference density stays fixed with cached geometry;
arbitrary variable density is always retained in the actual solve.
