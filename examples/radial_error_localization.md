# Fixed physical bands and cancellation correlations

Reduces already saved independent per-field-cell squared integrals and full
error cross integrals. Uses the same physical axis and wall bands, each 1/32
of radius, for every mesh; changing fluid element size cannot change the
reported physical region. Requires band faces to align exactly with cell
faces, and checks each squared-error partition sums to the original total.

Saves quadrature16/20 fractions for every actual term and error correlations
for the earlier full term audit, preserving its term order. Zero error-norm
correlations are explicitly placeholders, not measured correlations. No field
or derivative is created. Failed original cases retain their status and only
their driving/momentum evidence is reduced. Source and exact input hashes are
frozen; fresh replay verifies all metrics. This does not replace face jump,
local conservation or full trajectory accuracy controls.
