# Separately refined original-continuum rate study

Evaluates both signed physical perturbations through the complete wall-conforming
RHS on explicitly independent fluid/field meshes. One pair of actual rates is
reconstructed with independent quadrature16/20. Reports all nine physical-rate
components, full electric and current rate, initial representation errors,
axis/wall/max norms and complete induction/power diagnostics. No field, rate,
state, density or heat is corrected to match the reference.

Raw archives store complete initial/perturbed conservative states, magnetic
coefficients and both signed fluid/magnetic/current/electric rates. Expanded
physical samples are deterministically reconstructible; dtype/shape/value
fingerprints record them without duplicating hundreds of megabytes per case.
Every coefficient and physical-field fingerprint must replay bitwise. Stable
electric/current differencing occurs before their linear reconstruction, so
last-bit differences from the earlier sample-then-subtract audit are expected;
there is no change to the evolved or instantaneous physical equations.

`calculate(spec, epsilon, orders)` performs one complete check; `run` freezes a
caller-supplied explicit list and refuses overwrite. Spatial resolution changes
must be declared in the spec, including epsilon controls. Runtime is excluded
only when comparing replay metrics. These rates are not full trajectories.
