# Variational interface work

`variational_trace` evaluates the boundary part of
`K u_induced - K_fluid psi_total`, where `K_fluid` is the gradient of the
positive staggered fluid magnetic energy. This is the discrete natural
boundary load, conjugate to interface psi; no global energy residual is used.
At noncorner wall/end nodes it directly gives the dual-boundary integral of
tangential B. Interior residuals, including analytic-source discretization
error, are not folded into a boundary correction.

At an outer/end corner only the weighted sum of its two normal traces is
determined. Independent second-order exterior derivatives supply a consistent
initial split. Project them onto that single weak constraint with minimum
boundary-length-weighted squared trace change. Constant fields are unchanged;
smooth accuracy must be checked independently. This is a local boundary
Neumann-trace representation, not a thermal-energy repair.

`BoundaryPoynting.face_fields` supplies the corresponding Maxwell-stress
traces. `energy_fluxes` integrates the nodal product of **actual selected CT
electric field** and tangential B on each gas boundary, with annular dual
weights. The gas wall is impermeable, so it has no advective gas-energy flux.
Axis E and Br vanish regularly. Physical and numerical electric terms use the
same work flux; neither is counted as a separate added total-energy source.

Local fallback must retain the same boundary load for both reconstruction
orders and reevaluate its Poynting flux after corner selection. A troubled
cell's boundary face depends only on its own two corner electric fields,
preserving the checked-low-update contract.

## Inner wall of an annular chamber
For an inner wall (`mapping.r_offset > 0`), `variational_trace` also returns
`bz_inner = load[0]/dz`. The inside vacuum edge contributes +dz Bz to the wall
node's load: its outward normal is −r, and the sign is the opposite of the outer
wall. The two inner corners use the same minimum-norm corner split as the outer
ones, with one-sided stencils into the inner vacuum.

`face_fields` returns four entries. `energy_fluxes` returns a fourth, the +r
Poynting flux through the inner faces, E_θ B_z,inner/μ0. FluxUpdate puts it on
radial face 0 and subtracts it from the outward energy. The energy balance
closes through the inner wall: 6e-10 relative in
`examples/annulus_verification.py`.
