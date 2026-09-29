# Finite-time physical entropy verification protocol

Nine full-axis-transit shear-transport trajectories: nr8/16/32, nr16 time-half,
quadrature12, reversed current, larger entropy amplitude, ion gamma, and degree3.
The current is obtained from discrete Ampere on an independently initialized
physical field; no analytic current is substituted in evolution. Actual positive
conservative density is held fixed explicitly for this transport fixture.

RK4 uses 5*nz steps per transit (twice for the time control). Every accepted step
records physical entropy, proxy, internal energy and min/max pressure; five full
snapshots plus the before-final state are retained for replay. No energy or
entropy reset occurs. Independent physical norms at quadrature12/16 distinguish
entropy/temperature-wave errors from full-pressure error and initial projection
from final drift. Exact Bessel/radial integrals independently measure conservation.

Fresh output/source archive only. All exceptions and incomplete work are
preserved and make the driver fail. The source fingerprint covers this driver,
production and reference/audit code. A completed transport is not a self-consistent
magnetic or bulk-plasma trajectory and does not relax the full moving requirement.
