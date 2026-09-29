# Read-only compatible wave audit

Checks source/archive integrity, exact stored reference inputs, conservative
velocity reconstruction, Ampere current, affine bias and rounded export, all
reported physical norms and component energies. Original full Ohm, rate Ampere
and Faraday equations are replayed separately from the solver's residuals.
Physical norm/energy interpolation uses independent root-product polynomials.
The actual projected numerical bulk/Hall forcing and analytic pressure potential
are replayed from the stored physical inputs. Affine increment magnetic energy
is separately integrated and normalized by that small energy itself; the large
background is not permitted to hide a wrong weak-field energy.
Integration differences are measured rather than silently identified with
assembly quadrature. Full-field rate-relative conditioning stays explicit.

Rounded-attempt source is checked against its own immutable snapshot; only the
documented measurement/norm representation edits may differ from the current
source. Its fields are not regenerated or repaired. No trajectory, nonlinear
pressure or formation gate is inferred from these instantaneous audits.
