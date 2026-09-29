# Compatible physical wave-rate measurements

25 instantaneous linear cases: all three branches over nr8/16/32/64 at degree5,
four cubic slow controls, four exactly matched perturbation-only weak linear
controls, two physical quadrature controls, weaker full bias, original n1e20,
and opposite propagation. Main cases retain the full uniform bias in the saved
magnetic state as an explicit affine harmonic scalar and physical masses.
Both increment cochains and the lossy rounded-total export are saved. The exact
Bessel root boundary sets wall increments to zero before storage; no numerical
current is used to fit that boundary or velocity. Older case configurations
without an affine representation retain their original rounded-total semantics.
Exact parent fields supply only initial
cell averages/nodes and comparison targets, never the numerical bulk/current
forcing. Current follows Ampere from the new physical Hodge matrices.

This is the exact constant-background linearization, not nonlinear evolution:
momentum uses rho0 v, B0 in Ohm products and uniform inertia. A separately
represented pressure potential restores the complete linear E and original Ohm
and Faraday backward/rate-relative checks. Every field, forcing and input is
saved. Independent root-product interpolation at quadrature orders 10/14 checks
physical B/J/Bdot/Jdot/E, axis/max errors and actual magnetic/current energies.
Actual representation degree and integration order are recorded separately.

Original full Ampere backward error is checked against the un-cancelled physical
operator terms. Its error relative to the current after lossy total-field export
is separately retained; the affine representation does not pretend that export
can preserve every weak signal. Global harmonic/increment energy components
and independently integrated full local energy must agree.

The driver freezes source before any run and rejects directory reuse. An
unexpected exception remains a visible terminal error with prior records intact;
it is never mislabeled as physical initialization failure. No trajectory,
nonlinear energy/thermal coupling, variable density or formation claim follows.
