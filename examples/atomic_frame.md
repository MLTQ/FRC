# Heated argon in fixed reference coordinates

Continues the concrete failure in the boosted 3.5 µs heated argon source.
Two complete 1400-step histories use the same physical velocities, density,
temperatures, 20 eV/nucleus heat waveform and atomic inputs as that failed run.
One fixed frame equals the common physical boost [20000,-10000,5000] m/s;
the other differs by [137,-211,93] m/s to exercise a distinct state representation.
Neither the frame nor the physical conditions adapt during evolution.

The existing completed unboosted heated history is the comparison. Scalar
populations and temperatures must be frame-independent within numerical error.
Expected species lab momenta/energies follow the explicit Galilean transform
of those archived states. Expected current includes the actual charge residual.
This is coordinate consistency of the physical closure, not independent evidence
for its atomic rates or large-drift accuracy. Independent DI/RR/decay histories
and time-halving results remain in the original atomic_moments examples.

The existing `evolve` supplies full histories and independent conservation
measurements; no new test suite or generic reporting machinery is introduced.
Outputs save both the integration-frame and reconstructed lab states, fixed U,
private input paths/hashes, frozen source and comparison figures. No raw provider
tables or rate matrices are bundled. Failed prefixes remain explicitly failed.
Spatial/current coupling, photon momentum/absorption and atomic-data gaps remain.
