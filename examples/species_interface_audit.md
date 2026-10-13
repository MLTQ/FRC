# Complete independent weak species-entropy budget

## Purpose
Audit the entire current-driven linear entropy coefficient vector about the
archived stationary base, for both the saved-error and smooth independent mode.
It compares original production K with independent volume integration and the
full integration-by-parts decomposition in `species_weak_sampling`.

## Protocol
Orders 12, 16 and 20 distinguish original quadrature from refined independent
quadrature. Each order retains both the production q12 projection and its own
independently integrated L2 projection; neither replaces a production state.
The unchanged `RegularSpecies.transport` provides the original comparison.
Its fields, gradient and transpose equal those of the original wall wrapper;
the wall wrapper changes magnetic curl, which this diagnostic never evaluates.

All components receive the original -a/(e V) species multiplier and physical
entropy-rate reconstruction S(Cdot)/(S n0). Record component norms, signed
products, projections on the physical direction, and whole-vector discrepancies.
All three archived centered amplitudes are compared both in coefficients and
through the actual perturbed physical quotient. Both signs remain inputs.

The difference S(C/n) versus (S C)/(S n), and both radial/axial derivatives,
is measured explicitly. Strong advection from either representation is only a
comparison; the complete weak action retains one-sided faces and transfer.
Independent one-sided traces and volume-current fingerprints reproduce the
completed q16/20 charge audit. The extra q12 bridge is separately identified.

An analytic control uses Jr=r*k*cos(kz), Jz=-2*sin(kz), divJ=0 and entropy
r²/R²+0.2*cos(kz). Values, derivatives and the strong advective target are
analytic. It verifies the full-vector integration-by-parts quadrature and is
not subject to the physical fixture's no-normal-current wall condition.

## Contracts
Only the stationary current-driven entropy block is audited. No induction or
time solve, interface replacement, filter, floor, diffusion or mass change.
Native and physical norms remain distinct. The declared 1e-9 identity reporting
threshold concerns this diagnostic identity only, not a relaxed solver bound.
Raw vectors, projection matrices, all metrics/input hashes and source closure
are saved. Run through `guarded_species_interface`, then a full fresh replay.
No unit tests. Freeze calculation source after its first execution.
