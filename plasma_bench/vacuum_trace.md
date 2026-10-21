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

`interface_resistivity(..., speed_cap=None)` optionally bounds the local signal
speed in the LLF-scale interface resistance.
- **Why:** in near-vacuum wall cells (density-closure plasma) the fast speed
  reaches ~1e6 m/s. The resulting numerical dissipation then exceeded the coil
  work in NUCTE-III formation runs.
- **Effect with a cap of 2.4e5 m/s** (the fill Alfvén speed): the formation
  phase is unchanged, and the numerical heat falls from 85% to 29% of the work
  at 5 µs.
- The default is unchanged.

`tangential_field(mapping, flux, source=None)`: with the analytic coil `Field`,
`flux` must be the *induced* flux only. Its smooth harmonic exterior is
differenced and the coil B is added analytically: bz at the wall nodes, and br
at the end rows on the radial centroids.
- **Why:** in NUCTE-III case 6 the end-strap filaments (r = 0.1267 m) lie inside
  the two-interval wall stencil (wall + 2dr = 0.129 m at 16×128).
- **Effect:** the 1 µs response residual fell from 2.9e-2 to 1.0e-3.
- **What it did not fix:** the later growth of the case-6 residual, which
  converges with resolution (see `docs/NUCTE_BENCHMARK.md`).

## Inner wall of an annular chamber
`tangential_field` returns a 4-tuple `(bz_wall, br_bottom, br_top, bz_inner)`.
`bz_inner` is the vacuum-side Bz on the inner wall faces, from the one-sided
second-order derivative toward the axis, (3ψ_i0 − 4ψ_i0−1 + ψ_i0−2)/(2 dr r_in).
It is None without an inner wall. With an analytic source, the coil Bz at the
inner wall nodes is added, as at the outer wall.

`interface_resistivity` also marks the inner wall row. End-row marks start at
row max(i0, 1).

2.5D: `interface_resistivity` includes B_θ = √⟨r²⟩ g in the cell signal
speed. No tangential B_θ trace is needed: the vacuum B_θ is known exactly,
r B_θ = 0, or the prescribed threading value (`vacuum_flux.md`).
