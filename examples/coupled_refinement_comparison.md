# Complete physical comparison across F8/F16/F32

Requires both fine coupled probes and audits to have full fresh replays, plus
the original F8 probe/audit replays. Reconstructs all actual centered physical
responses independently on q16/20 and reproduces every archived response
fingerprint. All three grid pairs use differences of actual fields, with raw
squared-difference, reference, response and signed-product integrals.

Also reconstructs actual signed input currents using the original weak Ampere
map. Fine current coefficients must match the archived physical audits exactly.
Instantaneous physical conservation is recomputed for all three grids and
every sign/amplitude; every previously archived fine power metric must match.
This supplies the missing F8 physical-power baseline without a new RHS call.

No coordinate projection is treated as a common cross-grid eigenvalue. All
physical components, both species, finite-amplitude and tiny-signal effects are
retained. Native grid powers are separate from physical integrals. An
instantaneous or global comparison does not prove local flux/trajectory accuracy.
Save complete hashes/raw integrals and fresh replay under the original guard.
No solver, time evolution, changes to old evidence, or unit tests.
