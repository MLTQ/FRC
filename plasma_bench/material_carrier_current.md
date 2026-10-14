# Finite-mass carrier current on moving material

A coupled invariant sector: axisymmetric meridional fluid/current, toroidal B,
no bulk rotation, constant separate species entropies, zero resistivity.
All finite-mass Hall/self-drift, changing-density inertia, vorticity/stretching,
magnetic force and bulk drift stress remain. Constant entropy makes species
relative entropy transport zero; barotropic pressure enters induction only as
a curl-free potential. Equal temperatures eliminate Coulomb exchange. This is
not general thermal/coil/formation physics, and full electric scalar potentials
are not returned. The reduced electric field has the correct curl.

Physical pullbacks are G=F^-T G0 and B=F B0/J. Fixed material dN gives
Mn=integral dN G^T G. Let H=integral Bbasis^T Bbasis dV/mu0 and
D=-integral Bbasis^T curl(G) dV/mu0. Ampere is e Mn w=-D^T b, Jphysical=en Gw.
Thus current vanishes with carriers without dividing a reconstructed current
by a density floor. No exterior current/density/mass exists. The compatible
nodal gradient is in the electric space; its curl vanishes and Ampere has zero
weak charge work. Pointwise current continuity is not asserted.

Actual moving Faraday is H bdot=D E*. Eliminating finite inertia gives
(H+mu/e² D Mn^-1 D^T) bdot = D Mn^-1 rhs,
rhs=f - mu/e² Ddot^T b - mu/e Mndot w,
f=integral dN G^T [delta w x (B+mu/e curl w)-mu/e w x curl v].
Here mu=mi*me/(mi+me), delta=(mi-me)/(mi+me). Physical masses are retained.
The current equation retains both geometry derivatives, with no frozen metric.
A positive diagonally scaled direct solve enforces the actual equations.

Magnetic force is minus the fixed-flux energy geometry derivative. Drift stress
is +integral mu dN ww:grad(delta x). Their work cancels the exact geometric
magnetic/current exchanges. Species pressure and consistent fluid inertia come
from MaterialTwoSpecies. RK4 advances the simultaneous fluid and magnetic state.
No force/heat/energy repair, density floor, numerical diffusion or damping.

The initial/continuing entropy and zero-rotation restrictions are explicit.
Sharp thermal/kinetic edge validity, full species transport, nonzero resistance,
poloidal fields and actual finite-coil/vacuum coupling remain original-goal work.

The first fine free run stopped at step554 when electric-field reconstruction
lost Faraday accuracy (1.031e-10 against unchanged1e-10 limit). The source and
553 accepted states are preserved. The solver now refines the original mixed
rows H*bdot-D*E=0 and (mu/e²)D^T*bdot+Mn*E=rhs using the same Schur operator.
Corrections apply to the linear unknowns only, up to four passes. The original
full Faraday/Ohm acceptance checks remain unchanged; no state, force or energy
correction is made. Complete repeated fine trajectories finish 2000 steps;
maximum Faraday residual is below 3.04e-13 at the unchanged 1e-10 bound.

`project_initial` minimizes the actual magnetic plus carrier-velocity error for
supplied mutually consistent analytic B,J data. Its positive normal matrix is
H+mu/e² D Mn^-1 D^T. The right side is the magnetic target load minus
(mu/e)D Mn^-1 integral(dV G^T J_exact/e). No target J/n division is needed.
This constructs initial coefficients only; it never fits an evolved state or
forces its energy to a chosen total. Both B and J errors must still be measured.
The first free runs used magnetic-only projection and had large initial current
energy errors (85%/29% at4/8 and8/16); those runs are preserved. The new projection
reduces them to +3.49%/-0.152%, but increases B error and leaves 16.2% local J
error on8/16. Final free trajectories remain spatially unconverged. Read
`../docs/MATERIAL_CURRENT_COUPLING.md` for measured accuracy and validity limits.

MaterialThermalCurrent extends the sector with actual entropy/pressure work and
heat/exchange. This isentropic base and its active higher-order run stay unchanged.
Read MATERIAL_THERMAL_COUPLING.md for the thermal trajectory evidence and the
explicit upwind numerical diffusion introduced after a fine skew-entropy failure.
