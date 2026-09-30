# Shared particle continuity and physical flux

`ParticleFaces` holds transverse face averages qx=r n vr and qz=n vz, with
distinct periodic axial faces and explicit zero axis/wall radial faces.
`CompatibleParticleFlux.from_momentum` constructs these from the actual
conservative cell averages of n vr and n vz using degree-3/5 moment fits.
It is not a reconstruction from fitted numerical current or analytic velocity.
Normal faces use a single centered face stencil, addressed by integer face index.
The initial attempt reused a cell-interior two-stencil average and produced
positive periodic advection growth; its failed full-transit evidence is retained.
The corrected trace is distinct from the unchanged two-stencil DENSITY map.

`density_rate` is the sole finite-volume continuity law. `sample` lifts its face
flux to physical Gamma=(qx/r,0,qz) with a primitive in the normal direction and
the existing conservative density polynomial transversely. Consequently
div Gamma=-R n_t in each cell, with continuous normal flux and impermeable wall.
The quadrature points avoid r=0; regularity follows from qx(0)=0, not a floor.

`density_adjoint` and `vector_adjoint` return physical face work coefficients,
including all volume/metric weights. Fixed radial boundary variations are zero.
`integrated` converts to particle/s for independent cell/global budgets. These
adjoints are ingredients for momentum work; they do not themselves advance or
claim conservative momentum. No new density state, atmosphere, limiter or
numerical viscosity is introduced. All maps use the same immutable density map
as current and thermal energies. Radial physical Gauss order must be at least
2*degree+1: integration by parts pairs x-degree 2p products, hence r-degree
4p+1 after the cylindrical Jacobian. Underresolved quadrature rejects instead
of silently claiming entropy cancellation. This is a smooth positive-density
component.
