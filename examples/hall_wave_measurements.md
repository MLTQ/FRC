# Actual cylindrical Hall measurements

Runs the full declared case protocol and archives initial/final mode coefficients,
physical sampled fields, coordinates/weights and numerical radial eigenvalues.
Every accepted step records all three energies, positive heat, equation backward
error and a complex positive-k magnetic amplitude projection. Phase and damping
errors are measured from that evolved amplitude, not from eigenvalues alone.

Independent continuum Bessel fields supply physical-volume field/flow/weak-current
errors, true time-integrated heat and exact frequency. Independent higher-order
physical quadrature checks the modal energy metric and integration uncertainty.
Initial interpolation errors are reported alongside final errors. Time controls
can contain unresolved secondary numerical modes, especially on the slow branch;
a small energy residual is not a resolved-frequency claim.

Full configuration, dependencies, relevant production/reference/measurement source
hashes and lossless NPZ hashes are stored. The source fingerprint is restricted
to actual dependencies so unrelated future modules do not invalidate this record.
The production source must remain unchanged during a run. No reuse cache hides
reruns; reporting validates every case and file separately.

A separate helicity-space dense exponential retains every initial spatial mode
and supplies a semidiscrete temporal error in the complete energy norm. Three
extra slow-branch refinements resolve its previously exposed fast secondary
components. The oracle source and endpoint coefficients are also archived.
The residual norm includes all complex roundoff; it does not demand relative
real-field conjugacy of a near-zero difference of two accepted states.
