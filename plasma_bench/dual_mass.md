# Conservative primal/dual mass overlap

`DualMass` shares `DualCurrent`'s full-domain dual volumes in x=r²/2,z.
`density` maps a primal density (or signed rate) by four quarter-cell mass
contributions. `fluxes` integrates the lowest-order rectangular H(div)/RT0
interpolant matching the supplied **integrated** primal normal face fluxes.
The map commutes with finite-volume divergence, including the axis, radial
wall and periodic seam. It does not reconstruct a new velocity or suppress
the fluid method's numerical mass diffusion.

`primal_average` is the exact average of a nodal Q1 scalar. It is the volume
adjoint of `density`: sum(Vp q P f) = sum(Vd Q q f). This gives one explicit
electron energy convention. It does not assert that overlap densities are
exact point samples. `centers_x` gives actual dual-volume centroids for bulk
transport reconstruction; boundary centroids differ from the scalar nodes.

Radial input is (nr+1,nz); axial is (nr,nz+1), with **identical** periodic
endpoints. Output radial is (nr+2,nz), axial (nr+1,nz); dual axial index j is
the upper face of dual cell j. Flux units are preserved. Any finite signed
density/rate is accepted, leaving positivity to the species state. General
boundary fluxes are mapped; impermeability is enforced by the fluid caller.

Independent affine-face integration and nonuniform local divergence tests
are in `test_dual_mass.py`. Thermodynamic use is in `two_temperature.py`.
