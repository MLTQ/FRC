# Toroidal induction and scalar resistive work

Poloidal UCT uses full three-component fast speeds and the established staggered
magnetic jumps. It preserves the full one-dimensional LLF limit, axis gauge,
periodic axial endpoint and fixed conducting-wall flux.

`toroidal_resistivity` uses dg/dt=diff(Ez)/annular_measure-d_z(Er/r),
Er/r=-eta/mu0 d_z g. A compatible tridiagonal current mass P solves
P Ez=(eta/mu0) difference(s g), s=<r²>, with the regular axis equation and
Ez(R)=0 conducting wall. Ez(axis) is solved and may be nonzero. The same P
is the work metric of recovered Jz. Its positive node/edge-square decomposition
makes heat nonnegative; the semidiscrete toroidal magnetic-energy loss matches
that heat without an energy reset. P and its bounded geometry cache are shared
by the solve and diagnostic. Polynomial g=a+b r² with its conducting trace
has exactly the expected 8Db diffusion at every cell, including the axis.

Radial/axial Poynting fluxes use -Ez Btheta/mu0 and Er Btheta/mu0. Total energy
receives their divergence only; physical heat is diagnosed separately to avoid
double counting. Toroidal flux can change through the finite axis Ez even with
a conducting outer wall. The axis is not a second conducting material wall.
The discrete current weights are consistent numerical quadrature, not a claim
that every weight is an exact physical dual-cell volume. See the derivation and
independent continuum refinement in TOROIDAL_MHD_NUMERICS / VERIFICATION.
