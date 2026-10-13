# Physical spatial tangent and amplitude audit

`run(destination, probe_directory)` validates the completed probe, its source,
all inputs and full returned rate arrays. Rejected signs remain recorded and
receive no invented physical rate. Only the fluid node layout, field mesh,
independent root-polynomial/Fourier sampler and original analytic continuum
reference are constructed; there is no new plasma solve or time integrator.

`expand` differentiates the actual conservative-to-physical reconstruction:
density, all velocities, both specific entropies, B and both species pressures.
Pressure uses the original EOS and `dp=p*(gamma*dn/n+ds)`. Sampled density and
temperatures remain checked. `components` keeps scalar physical components
separate, including all current and electric directions. No unlike-unit norms
are added. `compare` reports global, axis-element, wall-element and maximum
differences; an exactly zero reference yields an undefined relative norm
(`null`) plus the actual absolute maximum, without a normalization floor.

The full physical-RHS derivative is the centered difference of physical rates
at the actual signed states. It is compared separately to the conservative
Jacobian action pushed through the physical reconstruction at the base. Their
difference includes finite-amplitude effects and the discrete base-stationarity
defect through the reconstruction's second derivative; they are not silently
identified. Current-rate differencing is also compared against the unchanged
weak Ampere map applied to the centered magnetic action. E is a solved algebraic
response, not an electric-field time derivative.

For each physical state component, report its directional projection rate,
orthogonal residual and residual against the separately declared coordinate
projection. These are directional diagnostics, not eigenvalues. Base physical
rates and their fingerprints remain explicit. The independent-mode control
retains original continuum rate, initial representation, both pressure,
current and electric comparisons. Every accepted amplitude's actual fields
are compared with its adjacent accepted amplitude control using physical differences,
not differences of error norms.

Quadrature16/20 records all component squared integrals/products, positivity,
full-field fingerprints and hashes. Raw coefficient inputs reconstruct every
field; full independent fresh replay must compare every metric and raw integral
array. Run through guarded_tangent_workflow only after collecting the physical
probe and its complete replay. No earlier source or state is changed, and no
diffusion, filtering, floors, artificial mass, atmosphere or energy correction
is introduced. Full trajectory, local budget and formation gates remain open.
