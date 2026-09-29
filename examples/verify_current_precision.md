# Increment precision verification protocol

Preserves the 41-case baseline and its rejected stationary controls. New records
carry the changed engine identity and extra audit/helper source hashes. Runs
18 repeated stationary controls plus 15 nonlinear trajectories: original pressure
and base wave meshes 8/16/32, wave amplitudes 1e-9 and 1e-12 T m on those meshes,
and weak strong-bias, reversed-helicity and strong-inertia controls at nr=16.

Reuses the independent physical energy/entropy, continuum wave and Faraday audit
on reloaded archives. Every case must finish all prescribed steps; weak-wave
field and current errors must refine relative to the perturbation itself.
Interrupted protocols reuse completed evolution records only after matching
the full configuration, engine provenance and archive hash; every audit is rerun.
Baseline metric changes include exact source and record hashes, keeping old
measurements distinct from reruns. Separate stored-state and nonlinear residuals
are retained. The 100-fold bias control tests storage limits without suggesting
that physically tiny waves are resolved below floating-point precision.
