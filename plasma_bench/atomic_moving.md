# Moving variable-species volumes and electrostatic vacuum

Current measured status: the 32-cell electron–ion run completes, but 64/128-cell
runs reject negative thermal energy in numerically diffused ion tails. Neutral
front position is not converged, and charged-front error remains large. This
research component is NOT adopted by the formation solver. See
`docs/ATOMIC_MOVING_FIELDS.md` for complete outcomes and the next physical fixes.

`MovingAtomicElectrostatic` stores moving axial faces x_j and integrated species
amounts Q_i=integral(cell state/Nref dz), with one fixed Nref. Each retained atomic
level and electrons keeps its own number, three momenta and total energy. Species
amounts advect relative to the mesh and optional local reactions can change them;
the mesh does not assume fixed electron/ion counts or two entropy variables.

## Geometry, transport and free boundary

Default interior mesh velocities interpolate between the two front speeds using
the actual face coordinate, preserving relative cell widths without remapping.
This is affine COORDINATE motion only; each physical species velocity remains an
independent unknown with full ALE flux. Optional `mesh_motion="barycentric"`
instead interpolates actual mass-weighted cell velocity. The latter allows outer
cells to stretch drastically when electrons outrun ions; the first complete
exploration exposed unresolved front fields and negative electron thermal energy.
The two outer faces move at the outward envelope of the occupied species' Euler
vacuum-fan speeds u±2c_s/(gamma-1), using the adjacent occupied cell state.
This first-order boundary Riemann state is explicit: a zero-density reconstructed
trace has no defined p/n sound speed and cannot supply its own fan. The exact vacuum trace has zero density,
pressure and ALE flux. No material cells exist outside those faces; `sample`
returns exact zero there. A conservative Riemann flux handles interior faces in
their rest frame (species_ale_flux.md). This tracks the OUTERMOST common support;
different species' internal fronts can still diffuse within that support.

Conservative number/momentum/energy reconstruction uses minmod gradients on the
nonuniform mesh and one-sided boundary gradients. One common slope factor makes
each species' traces admissible while retaining its actual cell mean; see
species_reconstruction.md. All reconstruction calls (including timestep queries)
count occupied and limited species. Evolved cell amounts are never clipped or
floored. HLLC/HLL/minmod diffusion remains numerical,
not physical transport. There is no mesh remapping; cell crossing or reaching
the fixed remote field ends rejects. Front accuracy must be measured separately
from exact emptiness outside the numerical front. This is 1D, not the existing
axisymmetric formation interface or a kinetic plasma edge.

## Field and work on moving volumes

Let A_s=Q_i[number_s], h_i=x_(i+1)-x_i, and let F_s denote the actual numerical
ALE number flux. The normal electric field at moving faces is

    E_(j+1)-E_j = e Nref/epsilon0 sum_s Z_s A_s,j.

E_left is a fixed input; E_right follows actual total charge. Between faces E is
linear, and it is constant in each empty exterior interval up to the fixed remote
ends a,b. Neither charge subtraction nor a neutralizing background is introduced.
The electric energy is integrated exactly for this represented field:

    U_E/(e Nref) = epsilon0/(6 e Nref) sum_i h_i(EL²+EL ER+ER²)
       + epsilon0/(2 e Nref)[EL_boundary²(x0-a)+ER_boundary²(b-xN)].

With mesh-face speeds wL,wR, continuity reconstructs the cell's physical number
current as linear interpolation of ALE flux plus (A_s/h_i) w(z). Integrating its
work with the linear E gives the cell energy source

    W_s = Z_s/6 {h_i[(2EL+ER)F_s,L+(EL+2ER)F_s,R]
                 + A_s[(2EL+ER)wL+(EL+2ER)wR]}.

Momentum gains e Z_s A_s(EL+ER)/(2 p0). Fixed-frame energy subtracts U dot force.
The sum of W cancels the independently differentiated interior PLUS exterior
field energy, including moving-volume terms. This is a semidiscrete continuity
identity; no global energy residual changes any state. It contains numerical
transport contributions and is not solely physical heating. Constant normal
fields at the remote ends can exert external Maxwell stress; momentum is only
expected conserved when their squared magnitudes agree.

## Integration and scope

RK4 evolves geometry and integrated state together. The step resolves relative
acoustic CFL, mesh contraction and physical local plasma frequency. Actual
stage thermal energies/Jacobians validate; rejection is retained without repair.
Optional local reaction/collision halves act on Q/h at fixed Nref. Photon/input
heat amounts are producer-cell bookkeeping; photon propagation is not implemented.
Sources and transport still require independent coupled time sensitivity.

This is the longitudinal moving species/field component. Transverse moving
Maxwell curls, magnetic forces, cylindrical geometry, physical coil/open-end/
wall conditions, nonlocal transport and efficient dense-plasma time treatment
remain necessary. It must not be described as a complete formation solver.
