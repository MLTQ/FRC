# Array assembly of atomic moment transfer

`VectorizedAtomicMoments` inherits the physical state, fixed Galilean frame,
field recovery, exponential-midpoint integrator and observations from
`FixedFrameAtomicMoments`. Only generator assembly changes: arrays of all
nonzero transitions replace the Python loop and its many small identity/index
allocations. Repeated origin/destination indices accumulate with `np.add.at`.
No reaction, state, density cut, population or exponential-series term is omitted.

The number/thermal/photon network block comes from the unchanged `AtomicNetwork`.
For ionization and equal-charge transfer, donor momenta and energies move with
the existing product/parent mass ratio; ionization supplies the physical
electron mass fraction. Recombination retains the existing centre-of-mass and
capture-energy equations in atomic_moments.md. Signed momentum and cross terms
remain signed; no positivity clipping, tolerance or budget correction is added.
Roundoff ordering can differ from scalar assembly. Existing missing channels,
rate validity, small-drift closure and photon recoil assumptions remain.

This is an explicit alternative source accepted by MaterialReactionOverlap and
ReactiveMaterial through their existing source contract. It does not change
spatial or temporal resolution, install coarsening, or make a full formation
claim. Compare complete trajectories and independent conversion solutions
before adoption; a faster generator alone does not resolve mesh proliferation.

Measured in outputs/atomic_vector_moments: all ten complete histories succeed.
The two full argon pairs take2.52x/2.33x less elapsed time with identical
populations and maximum current difference3.41e-12 A/m². Independent DI/RR
population and donor-moment errors decrease about4x on timestep halving.
This source is used explicitly in examples/reactive_material_fine; the original
source remains available. See the measured moving-field status for full-run limits.
