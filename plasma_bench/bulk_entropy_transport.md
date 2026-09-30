# Prescribed bulk flow with evolving density and two entropies

`BulkEntropyState` stores actual number-density cell averages and two nodal
specific log entropies. `PrescribedBulkTransport` advances all three together
with the same conservative face continuity and physical scalar metric. Its
callback supplies cell means of particle momentum from the current numerical
density and an explicitly prescribed velocity. Momentum, current and magnetic
field are NOT evolved; this fixture is not a force-balanced moving plasma.

Both species pressure states validate at every RK4 stage and the endpoint.
Failure leaves the caller's arrays intact, with no retry, floor or energy reset.
`thermal_work` integrates the independently exposed face pressure-work rate with
the same RK4 weights. Comparing this ledger with physical internal-energy change
measures finite-step error; neither quantity is set equal to the other.

Physical entropy is a bilinear function of evolving density and nodal entropy.
RK4 does not exactly preserve general quadratic invariants, so finite-time
entropy drift must be measured and refined. Particle conservation is linear and
telescopes under the shared faces. This driver is explicitly a smooth bounded
verification fixture; it supplies no mechanical energy or complete momentum law.
