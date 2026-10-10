# Physical current and electric amplitude differences

Reads two terminal successful signed rate workflows on the same physical grid
and independent reference. Validates their complete source closures, reference
and raw-array hashes, and original declared induction bounds. It reconstructs
the centered current rate and E separately for each amplitude at every original
quadrature 16/20 point, using the independent polynomial sampler. Reconstructed
field, reference and weight fingerprints must match each original full audit.
Each response is `(plus - minus)/(2*epsilon)`, so absolute differences refer
to a unit-amplitude response, not to the actual small signed perturbation.

The volume norm of the physical difference is divided by the independent
continuum norm; normalization by each reconstructed field is also retained.
Current components and the vector norm are reported separately. Per-radial-row
squared integrals, maxima and quadrature differences remain available. These
are physical field differences, unlike changes of physical error norms or
Euclidean coefficient norms. They do not establish asymptotic amplitude order.

Only the original fluid sampling geometry and field mesh are constructed;
there is no field-mass factorization, plasma solve, new rate or state change.
Two radial field cells are reconstructed at a time, retaining every axial
point and the full scalar weight matrix. A complete postprocessing replay
must reproduce all arrays, fingerprints and non-runtime metrics. This is
end-to-end saved-workflow verification, with no unit tests. It measures J_t
and E only; it does not certify fluid/magnetic difference norms, evolution,
formation, vacuum coupling or confinement. Prior source files stay frozen.
