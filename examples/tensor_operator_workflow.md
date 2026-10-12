# Actual operators observed during complete physical workflows

Runs both full signed F4/64, B16/64 pilots with the unchanged independent
radial1,N64 reference. Observers fingerprint every actual base and strong-wall
mass, Hessian, curl, basis, quadrature, boundary projection and coordinate array
immediately after normal construction. No constructor or solve is skipped.
The original and tensor fingerprints must agree exactly.

Each full physical calculation must reproduce its archived pilot, including
every saved input/rate/face array and non-runtime metric. The observations have
a separate explicit source archive; the numerical workflow's own source record
remains accompanied by this observation record. Complete observation replay
checks the captured operator records as well as both physical reruns.

`fingerprint` includes sparse format, shape, dtype, indices and values.
`operators` records the assembled geometry used by the actual full simulation.
`run` restores all temporary serial observers even on failure. No state or
operator is mutated by an observer, and no unit tests are created or run.

Measured operator identity is limited to the recorded coarse grid. Fine-grid
acceptance additionally requires complete original/new physical comparisons,
strict residual/work checks and a fresh full replay; source-level separability
does not by itself certify accuracy, formation or hardware conditions.
