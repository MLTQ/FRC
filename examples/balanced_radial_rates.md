# Fixed radial16 rates after balanced strict refinement

Separate opt-in copy of the frozen radial_mode_rates protocol, substituting
BalancedWallPlasma only. Same independent N128 reference, input sampling,
quadrature16/20, full physical signed RHS and initial/face/current diagnostics.
All strict bounds and physical equations remain unchanged. Original accepted
and failed records are preserved. This source receives its own frozen closure
and full replay; no old archive is overwritten or reinterpreted as passing.

Captured-state recovery and analytic cylindrical polynomial controls justify
using this variant for further bounded instantaneous consistency measurements.
They do not establish full response, stability, exact vacuum or formation.
