# Full signed spatial workflow with streamed audit

Same independent radial16,N128 reference, fixture, exact signed conservative
inputs, SplitInverseWallPlasma, physical masses, solver bounds and complete
RHS as split_inverse_radial_rates. Only postprocessing volume arrays are
batched through streamed_rate_audit and streamed_interface_audit. Every
quadrature16/20 point and original physical field/face diagnostic remains.

The driver saves complete exact inputs before each solve and full returned
rates before expensive postprocessing, with an explicit progress status.
Completed case summaries include all original diagnostics and elapsed time.
Strict solver failures remain failures with exact input and error; auditing
exceptions or process interruptions leave visibly incomplete progress, never
a passing terminal case. Later attempts use fresh destinations.

Source/runtime closure freezes before initialization and is rechecked at the
end. Repeated runs compare all coefficient/face arrays and non-runtime metrics
bitwise. Original-versus-streamed comparisons retain diagnostic rounding
differences separately. chunk_cells controls audit batching, not mesh or solver.

The prior large attempt stopped during full-array postprocessing when the
host reported no writable disk space. That preserved interruption motivates
this change; it is not evidence of a physical solver rejection. These are
instantaneous consistency workflows, not complete trajectories or formation.
