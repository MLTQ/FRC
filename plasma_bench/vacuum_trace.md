# Exterior magnetic traction at the gas boundary

`tangential_field` obtains the vacuum-side tangential B using second-order
one-sided normal derivatives of the total flux at the gas surface. Shared normal
B remains the fluid CT face flux. Exterior tangential B supplies Maxwell stress
and Poynting flux at the boundary, including force/work on interface currents.
Using only reflected interior tangential B loses that work even for ideal flow.

Two exterior grid intervals are required (already enforced by CouplingMesh).
Corner values are averaged along the boundary faces with the same regular-axis
centroid convention. This is a finite-resolution interface approximation and
requires separate stress/work refinement checks; it does not impose continuity
of tangential B across an unresolved current sheet or move the gas surface.
`interface_resistivity` assigns positive numerical D=a h/4 at actual plasma
interface nodes. The half-cell dual area yields approximately E_num=a jump(Bt)/2,
the LLF boundary dissipation scale. Its sign uses the symmetric stiffness
current, avoiding indefinite one-sided derivative-jump diffusion. It is zero
at all interior plasma nodes, vacuum nodes and modes. Its electric field also
enters all adjacent Poynting fluxes, with separately diagnosed positive numerical
Joule work. This is numerical interface diffusion, not physical conductivity.
