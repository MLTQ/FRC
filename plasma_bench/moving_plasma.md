# Experimental material plasma with a genuine vacuum

`MovingPlasma` is an ideal, adiabatic, axisymmetric material-mesh reference.
Its unknowns are material vertex positions/velocities and fixed flux labels.
Cell mass and p V^gamma are fixed; nodal inertial masses are positive integrals
of initial physical density times Q1 shape functions. No exterior fluid arrays,
density/pressure floors, artificial viscosity or thermodynamic repairs exist.
`MovingState` copies coordinate and velocity arrays independently.

The potential is gas internal energy plus full-space induced magnetic energy.
`evaluate` combines exact discrete volume derivatives, magnetic geometric
virtual work and the material-flux chain rule. Vacuum coefficients satisfy
their stationary elliptic equations at each force evaluation. A prescribed
uniform axial field supplies `-0.5 Bdot sum(material magnetic load)` source
power. Scalar bias is in tesla; optional waveform objects have amplitudes
interpreted as tesla here, with an explicit source time cap required.

`step` uses velocity Verlet, second order in time; all proposed geometry/field
and velocity values validate before state, clock or work is committed. The
axis radial reaction is constrained and performs no work. Energy is diagnosed,
not forced exact by a correction. `advance` uses an explicit material-wave and
relative-node-motion cap plus any declared source cap. Its optional callback
observes accepted steps only. `diagnostics` distinguishes response work, true
physical mass, mesh quality, boundary motion and vacuum residual.

This reference currently has no resistivity/heating, ionization/transport,
arbitrary spatial coil sources, remapping, shock capture or open gas ends.
Q1/P0 material discretization can have low-energy mesh modes and distortion;
positivity of quadrature Jacobians is not a global topology/robustness theorem.
Uniform initial density used with the analytic FRC pressure is a mathematical
verification choice, not a physical species/temperature closure. Independent
force, expansion, equilibrium, mesh/time and source-work checks are required
before interpreting evolution. The existing Stage 2 formation gate is unchanged.
