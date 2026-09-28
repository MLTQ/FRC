# Quantified verification report
`vacuum_metrics` measures field reconstruction and face divergence. Manufactured
Bessel/cosine modes separate diffusion spatial and temporal errors.
`experiment_sensitivity` varies time step, mesh and exterior-domain size at
otherwise fixed RLC conditions, with topology disabled for performance and
explicitly labeled. `main` writes numeric metrics and a human-readable report.
All errors use explicit references; no result is fitted to produce reversal.
Two additional refined oscillatory runs test sensitivity of the topology flags
and sampled lifetime, keeping the snapshot interval fixed at 1 µs.
