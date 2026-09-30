# Cylindrical joint transport prototype

Gauss volume points evolve actual n, radial momentum, angular momentum L,
axial momentum and both physical entropy densities. Their physical integration
weights are 2 pi r dr dz. Radius zero is an exact boundary, never a positive
radius floor. Radial boundary states are evaluated from interpolated physical
energy variables; density/temperature traces must remain positive.

Radial pair flux uses area bar(r) and angular flux bar(r) Gamma bar(vtheta).
For Q_N the paired integrated centrifugal radial source is
sum_j Q_ij bar(r) n^log [r] bar(vtheta) bar(Omega). Its work cancels the
extra angular kinetic-work jump, while L remains purely conservative.
The pressure/geometry combination is evaluated as a pressure difference,
−sum_j Q_ij bar(r)[p]/M, avoiding cancellation of a constant pressure source.
These are derived geometric terms, not a force fit or global energy correction.

Volume and entropy-projected boundary residuals are lifted with [I,Rᵀ]. A
shared centered interface flux conserves mass, axial/angular momentum and both
entropy sums. Axis/wall have exactly zero mass flux and a pressure-only normal
momentum flux. Periodic axial transport uses the existing high-order two-state
composition. The physical-energy spatial proof uses the actual angular H
gradient M Omega and the same node/trace states, not an independent account.

No artificial viscosity, heat, floor, post-step reset or atmosphere is added.
RK4 checks every stage and final boundary trace and rejects atomically. This
is a smooth positive-density experimental component, requiring original full
two-dimensional, near-axis, continuous-quadrature and time verification before
adoption. Gauss quadrature conservation is not exact continuous-energy accuracy;
axis regularity between nodes and boundary accuracy must be measured explicitly.
