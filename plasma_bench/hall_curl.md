# Hall electric quadrature and adjoint curl

Geometry-only experimental work pairing on the existing ToroidalState magnetic
variables. C maps electric coefficients to psidot/gdot with exact CT: fixed
axis/wall psi, Ez(R)=0 and periodic z. Electric bases are Er=r times constants
in each annulus and continuous linear z; Etheta continuous linear r/z with
zero radial endpoint traces; Ez linear in r² within each annulus and constant
in each z cell, zero at R. Ez(axis) is free. Three component interpolation
matrices I use common physical tensor-Gauss points and 2 pi r dr dz weights W.
The component masses M=I.T W I are positive; sparse factors are geometry-only.

`magnetic_gradient` differentiates exactly the existing cell-poloidal plus
annular <r²> g² energy. It does not change pressure recovery or repair energy.
The electric-space current solves M j=-C.T grad(H), so its work against any
electric coefficient vector is exactly minus the magnetic energy rate.
`current` samples that Riesz representative at the common quadrature points.
`project_electric` uses the same I/W/M, essential for Hall work cancellation.

`magnetic_samples` uses linear Br/Bz face interpolation and regular r g.
The production magnetic energy is still its documented cell-average convention,
not the squared norm of these point samples. Current consistency, local work
and boundary convergence must therefore be measured; algebraic cancellation
alone does not validate this stencil. The dual current inherits electric-space
traces, which can constrain boundary accuracy. This is an experimental spatial
component, not enabled in ToroidalMHD, and has no time integrator or fluid force.
`hall_poynting` separately audits the physical energy flux; it does not repair
the magnetic gradient or the local Hall transfer error.
