# Independent compatible-energy verification

Checks exact uniform-field energy and second-order convergence to an analytic
harmonic-polynomial magnetic-energy integral. A compact field compares the
fluid quadrature with the independently assembled exterior stiffness; a
nontrivial symmetric perturbation checks its variational gradient.

Independent harmonic-polynomial Br/Bz verify natural boundary traces, including
the corner projection and regular axis. The trace refinement includes unequal
radial/axial spacings, exercising distinct half-edge weights, not only square
cells. A prescribed electric curl in uniform current-free B verifies that
Poynting transport changes magnetic energy, not
thermal energy, cell by cell. These checks are distinct from simply observing
a small global response residual.

A short-wave flux fixture exposes energy missed by squaring averaged fields.
The new convention must reject reinterpretation of that evolved legacy state
without modifying it; explicitly constructing initial total energy from the
same gas state recovers the specified pressure. Full shock/source refinement
and time convergence are measured separately.

Actual first- and second-order fluxes at an evolved cylindrical shock must
satisfy the semidiscrete full-space response identity with uniform prescribed
bias. Independent symmetric differences of actual response energies after
vacuum recompletion check the same zero-power result without modifying state.

A manufactured boundary momentum deficit forces shared-corner fallback. Every
flagged cell retains its independently valid low-order derivative, and an
independent nodal dual-area sum verifies the actual mixed outward work and
integrated fluid energy loss. This exercises the new boundary-work path during
local fallback, rather than checking only all-high or all-low fluxes.
