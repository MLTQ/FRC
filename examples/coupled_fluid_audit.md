# Complete physical audit of the refined-fluid RHS

Separate extension of the frozen trajectory tangent audit. Requires the fine
probe's complete fresh replay and retains its full q16/20 physical RHS versus
conservative pushforward, both species pressures, independent continuum mode,
actual adjacent-amplitude differences, positivity, E and all J components.
Coordinates retain their fine-grid metric; physical component comparisons are
the authority for cross-grid interpretation, not coordinate projections.

Reconstructs the actual F8 signed RHS on the same physical grid and requires
all four baseline response fingerprints to reproduce the original completed
audit exactly. Reports actual fine-minus-F8 physical fields and their hashes.
No differences of error norms stand in for differences of fields.

Uses the unchanged weak Ampere map on each actual signed magnetic input to
archive its current coefficients; no floating homogeneous scaling is assumed.
`coupled_rate_physical.power` measures full physical instantaneous energy,
particle, entropy and momentum integral rates at the actual signed states and
base. Original grid powers remain separately reported. There is no correction,
local-flux acceptance or trajectory-energy claim.

The complete raw current table and physical integral table, scientific metrics,
source and inputs must replay exactly in a fresh process under the standard
resource guard. Old solver/audit sources are untouched. No unit tests.
