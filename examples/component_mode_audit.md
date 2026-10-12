# Independent component errors and physical localization

Reconstructs complete accepted signed rate workflows for two amplitudes on the
same grid/reference. Uses the original physical quotient derivative and root-
polynomial sampler. All reconstructed physical-rate, E, J_t, reference and
weight fingerprints must reproduce the original whole-workflow audits.

Reports all nine fluid/magnetic components, all three E and J_t components,
and separate full E/J_t vector norms. Each component uses its own independent
continuum norm. First/second errors and their actual physical difference are
integrated over the volume, r < 0.05 R, r >= 0.95 R and the remaining interior.
The fixed bands are diagnostic regions, not new physical boundaries. Squared
error fractions, row integrals and sampled maxima remain available. Zero
reference norms give an undefined relative norm rather than an invented scale.

No field-mass factorization or plasma solve occurs; all original quadrature
16/20 points remain, streamed by two radial field cells. This measures saved
full RHS accuracy and sensitivity, not time evolution or formation. Complete
reconstruction replay must match every array/fingerprint/non-runtime metric.
No unit tests, diffusion, filtering or repairs. Earlier sources remain frozen.
