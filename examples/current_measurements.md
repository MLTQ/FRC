# Complete current-inertia measurements

Runs the frozen 35-case protocol, recording complete inputs, every-step energy
and equation residuals, actual initial/final field/current/heat arrays and
SHA256/environment provenance. Reuse requires completed records with identical
configuration, code and NPZ hashes. Historical reports retain their old sources.

Independent continuum B/current/local-heat errors use physical volume weights.
Whole-space analytic energies use separate spherical quadrature, checked at
orders 48/96 and against the uniform decay's closed formula. Initial weak-current
error is retained. Source screening starts at zero and uses the exact prescribed
linear bias ramp; no numerical source circuit is assumed.

These records have no moving fluid or thermal feedback. Local carrier/thermal
validity and global energy consistency are separate requirements. No incomplete
case is reused or presented as a finished measurement.
