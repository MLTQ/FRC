# Paired sources for the moving plasma

`collision_config` validates the optional `maxwellian_pair` closure, explicit
positive Coulomb logarithm and monatomic species. `collision_sources` evaluates
one force/heat pair at actual positive plasma quadrature. No vacuum samples,
atmosphere, state floors, clipping or energy correction are introduced.

With solver drift w=Ui-Ue and J=e*n*w, the electron-relative drift passed to
`pair_moments` is -w. The returned electron force enters Ohm as R_e/(e*n), so
J.E_collision=R_e.w is the positive relative-flow energy lost to collisions.
The same call returns Qe and Qi, whose sum equals that loss. The ion force is
opposite; net bulk force is zero. Both entropy sources use the existing material
dual volumes, including negative species heat during temperature equilibration.

The default retains prescribed eta*J and the previous optional small-drift
thermal exchange, including the original mass-weighted Joule partition. Paired
mode requires eta=0 and separate exchange=False to avoid double counting.
It uses the exact instantaneous displaced-isotropic-Maxwellian moments described
in `maxwellian_pair_collisions.md`, not a conductivity tensor or kinetic solver.

The local relative relaxation rate is K/n*(1/me+1/mi), where R_e=K*w. It is a
local collision timescale, not a full-system timestep limit. The separate scalar
temperature-exchange frequency is null for paired mode; the actual exchange
power remains available. The explicitly supplied logarithm is reported as given,
not inferred from an edge floor or a claimed finite-drift screening model.

`MaterialVectorCurrent` uses the returned electric field and heat. Its driven
`MaterialCoilCurrent` child reuses that electric field in its electromagnetic
solve and inherits the already-applied heat, avoiding divergent closures or
duplicate heating. This isotropic collision law has no explicit B dependence.
Anisotropy, heat flux, nonlocality, magnetized corrections, a justified collision
cutoff and five-moment validity remain unresolved for formation trajectories.
