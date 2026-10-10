# Streamed full physical rate audit

Reconstructs the same full signed nonlinear rates and independent analytic
target as matched_mode_rates, at every original quadrature16/20 point. Two
radial field cells are processed at a time. The entire axial representation
and all nine physical components, E, current, initial representation and
axis/wall/max norms remain. This changes diagnostic memory use only.

NormAccumulator retains radial-row squared integrals and sums with math.fsum;
its reduction order differs from the older full-array np.sum. Matrix batching
can also change final rounding. Old/new differences must be measured by a
complete identical physical workflow, not assumed bitwise equal. Each new
workflow must separately replay bitwise. FieldFingerprints hashes each full
logical array's shape, dtype and ordered row bytes without retaining it.

measure uses the original independent TransportSampler, physical quotient
derivative and continuum target. The sampler still retains a full scalar
quadrature-weight matrix; only the much larger physical/gradient arrays are
streamed. No points, terms, modes or state coefficients are removed. No solver
equation, acceptance tolerance, diffusion or physical boundary changes.

Complete signed workflows and full replays provide verification under the
user's end-to-end-only preference. This is an instantaneous spatial audit,
not a time trajectory, formation calculation or confinement prediction.
