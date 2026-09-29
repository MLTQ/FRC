# Conservative dual current

`DualCurrent` maps primal g=Btheta/r to integrated poloidal charge currents
on scalar-node dual volumes. In x=r²/2, K=2pi <r²>g/mu0 is the current stream
function at primal cell centers. Radial current is K(z previous)-K(z current);
axial current is K(radial right)-K(radial left). The same differences cancel
local divergence, including the axis and periodic seam, without a projection.

The wall stream function uses the **explicit physical Btheta trace**. Omitting
it silently inserts a return-current sheet at the boundary. A z-dependent wall
trace yields nonzero local radial current, which cannot be discarded when
transporting electrons. Constant trace gives zero radial particle flux.

`gradient` uses scalar node values on x,z and supplies Er/r,Ez coefficients.
`faraday` is the physical commuting curl and retains Ez at both radial ends.
`work` verifies magnetic_rate = -current_work + boundary_magnetic_input.
The boundary input is -integral(S dot n); it must not be dropped when Ez at
the wall is nonzero. Restricted Ez(R)=0 recovers the existing Hall work pair.

Scalar values are nodal/finite-element values with lumped dual volumes, not
silently assumed exact dual cell averages. The total volume is the full cylinder;
the axis and wall each have their actual half-support volume. This component
does not solve fluid motion or select a physical thermal boundary model.
