# Simulation orchestration
`run(Experiment)` builds analytic sources and an independent optional stationary
solver, accumulates every-step scalar histories, and saves field/topology frames
at output_stride plus endpoints. Probes report nearest-node actual positions.
It separates vacuum counterfactual reversal, plasma-region reversal, on-axis
reversal (oriented by the bias sign) and closed
flux, with a 1e-6 relative-bias sign tolerance. Zero-bias ratios are NaN in numeric
archives but null in JSON; signed reversal and its phase classification are
undefined then. Topology may be
disabled for convergence benchmarks and is labeled unevaluated. Physical fields
use analytic external B plus reconstructed induced B. Energy histories refer to
the exact induced-field identity and separately to chamber-integrated total B.
Outputs explicitly refuse to establish FRC formation, stability or confinement.
