# Independent face-radial fluid checks

An analytic regular, wall-compatible polynomial force must be recovered in the
true variable-density radial metric, including the first annulus. Physical K
is checked against independent continuum quadrature. A fully 2-D radial-only
transport decomposition verifies the retained diagnostic total-work expression
with an independent higher-quadrature sampler. That diagnostic force is not
used by the evolved physical-force method; actual axial/angular momentum
fluxes remain conservative. The evolved method retains its convective physical
energy defect, separately from exact pressure/thermal work pairing.

Original 2-D Euler acceleration and curl references must converge with a
separate near-axis radial-force accuracy gate. A short unforced nonlinear wave
must refine and a deliberately rejected step must leave every state array
unchanged. Longer full-transit/equilibrium tests remain a separate protocol.

The actually evolved pressure force separately cancels physical thermal work;
that test makes no claim of exact convective total-energy conservation. An
interrupted nonlinear wave is serialized, hash checked and resumed; its entire
accepted history, state/rates, snapshots and energy ledger must equal an
uninterrupted integration bitwise, preserving the analytic-wave physical inputs.
