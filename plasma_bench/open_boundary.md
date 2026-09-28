# Variational infinite magnetic exterior

`multipoles` evaluates psi_l=(a/R)^l (1-x²) P_l'(x)/sqrt(l(l+1)), l>=1,
with R=hypot(r,z), x=z/R and analytic r/z derivatives. These regular-axis modes
solve Delta*psi=0 and decay at infinity. They represent only the induced field;
prescribed source coils/bias remain separate. The source-enclosing sphere must
lie inside the matching rectangle, enforced by the fluid coupling layer.

`exterior_energy` uses Green's identity to integrate the entire unbounded
exterior's Dirichlet energy on the three non-axis rectangle edges. Gauss
quadrature, reciprocity and a positive-eigenvalue check control assembly.

`MultipoleBoundary` couples those modes to the original interior energy stencil.
Unknowns are interior nodal flux followed by mode amplitudes. Boundary node
values are the modal trace. Interior/mode cross terms, boundary half-cell
tangential energies and the analytic exterior energy form one symmetric matrix.
The small dense modal block is symmetrized at roundoff after reciprocal assembly.
`pack` rejects boundary data outside the trace space; `expand` restores a full
nodal field. The matching surface is not a conducting wall or current sheet.
`outside_field` evaluates flux and analytic magnetic components beyond/on the
matching rectangle, including the regular axis limit. It rejects interior points
and mismatched amplitudes. No exterior field is extrapolated from a finite box.

Finite mode count is an approximation, and the nodal interpolation of the
analytic boundary trace introduces spatial consistency error. Both mesh and
mode refinement, and movement of the matching surface, must be verified.
Energy symmetry is not a convergence proof. The module adds no fluid outside
the chamber and makes no moving-free-interface or stability claim.
