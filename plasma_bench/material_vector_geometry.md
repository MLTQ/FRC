# Full material field/current geometry

`vector_geometry` applies the complete axisymmetric material deformation F,
including r*grad(theta) and r/r0, to electric one-forms G=F^-T*G0 and magnetic
two-forms B=F*B0/det(F). Curl(G) uses the same contravariant Piola map. No angle,
rotation, density or thermal state is replaced. The exact empty exterior belongs
to the vector space's compact-field boundary contract.

Mn=sum integral dN G_i.T*G_i, H=sum integral dV B_i.T*B_i/mu0, and
D=-sum integral dV B_i.T*curl(G)_i/mu0 use all three physical components.
Actual material geometry and both temperatures come from MaterialTwoSpecies.

In the rotating cylindrical frame, Fdot=L_eff*F with
L_eff=grad(v)-Omega*J_rotation. Its entries are meridional velocity derivatives,
r*grad(Omega), and vr/r. The omitted frame rotation is skew and performs no
metric work. Physical vorticity retains the cylindrical connection terms,
including (curl v)_z=2Omega+r*Omega_r; it must not use curl(L_eff) as a substitute.

The shared magnetic metric rate is L_eff+L_eff.T-div(v)*I. Hdot and Ddot both
use it. Electric-map rate is -L_eff.T*G; its contribution to Mndot is exact.
No frozen metric or missing angular geometry term is permitted. A full
current/fluid/torque trajectory must exercise these maps before a validity claim.
