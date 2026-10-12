# Resource-bounded component-stress signed calculation

Selects component_stress_rates inside the original guarded_tensor_refinement
driver, extending the explicit source identity and restoring serial bindings.
F64/256, fixed radial16/N128 reference, epsilon1e-4 and tolerance1e-10 remain.
Field radial/axial counts are explicit CLI inputs; --original requires a fresh
complete scientific replay of the same source/settings.

Keeps the 8 GiB own-process cumulative RSS and 256 MiB free-disk thresholds,
SIGINT and five-second self-exit fallback. Final values also need checking.
All partial records survive; resource interruption differs from strict numerical
rejection. No resource exception may weaken physical checks or certify accuracy.
No unit tests. This wrapper does not itself prove refinement feasibility.
