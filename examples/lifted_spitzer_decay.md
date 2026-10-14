# Same Spitzer experiment with magnetic coordinate addition

Reuses the complete physical setup, independent reference, observations and
unchanged acceptance bounds of spitzer_bessel_decay. The standalone entry
temporarily selects LiftedModalGauss in that module, restoring the original
binding afterward. The original experiment file and all archived failures
remain unchanged. This process-local selection is not intended for concurrent
threads inside one Python interpreter.

The provenance seed includes this driver, the original experiment and the new
integrator closure both before and after evolution. Output explicitly names
the integrator variant. No state, equation, coefficient, time step or reference
is changed by this driver. Runs require a fresh output directory and preserve
all accepted arrays/stages exactly as the original experiment does.

Purpose: complete the same full80/160-step physical evolutions after repairing
the boundary-coordinate cancellation mechanism, with original wall/solver
tolerances. Read plasma_bench/lifted_modal_gauss.md for the numerical change
and docs/SPITZER_DECAY_RESULTS.md for measured physical outcomes and limits.

The80/160-step full12us runs are complete, actual exits0 collected. Energy
drift decreases from4.65774e-4 to2.91367e-5 of initial field/current energy.
Magnetic and current errors do not improve with this time refinement; their
remaining spatial/nonlinear/numerical contributions are not claimed separated.
