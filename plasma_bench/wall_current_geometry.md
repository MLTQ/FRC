# Physical current at the driven wall

Adds one azimuthal nodal function per axial node. Write chi=phi_wall-P0 phi_wall,
where P0 is physical L2 projection into the old zero-wall theta space.
Thus chi(R)=1 and <chi,E0>=0. Current is J0(B)+chi*y(z); y is an evolved current,
not a fitted boundary value or diagnostic correction. The full current mass is
computed from the same physical quadrature. Magnetic geometry and admissible
Faraday curl remain unchanged; the extra E_ind coefficient is constrained zero.

Old coefficients embed with one additional theta row before Ez. `nodal_theta`
converts that orthogonal representation to actual physical nodal values for
independent reconstruction. Gradient potentials have zero theta component.
The new current differential uses the same polynomial field and its derivatives
on both sides of each interior face. All original advection/stress terms are
inherited without changes; their maps now include the wall mode.

Only filled periodic rigid-wall geometry is supported. This is not a vacuum,
plasma edge, sheath or finite-coil model. No current or heat filter is applied.
Old geometry, transport and archived simulation sources remain unchanged.

Measured full24µs160/320-step runs and their remaining spatial/current limits
are recorded in docs/WALL_CURRENT_RESULTS.md. The finer spatial run is still
active; this is not yet an accepted formation or vacuum-boundary model.
