# Adaptive physical diffusion substeps

`integrate_diffusion` is a pure fixed-geometry wrapper around midpoint resistive
steps. It compares one full step with two half steps in full-space magnetic
energy norm, physical-current weighted norm and local-heat-density norm. The
accepted solution is the uncorrected two-half-step solution and its sum of
positive physical heat/source work. Rejected trials change no accepted plasma
state. There is no extrapolation, pressure floor or energy repair.

Reference norm scales retain the largest accepted amplitude over the requested
span, so decayed roundoff tails do not require resolving forever. The initial
energy per square-root volume, using the kernel's SI permeability, supplies the
absolute heat-norm scale. These are
error normalizations, not field/density/pressure floors. Tolerance is a local
relative estimate, not a theorem of global error or magnetic monotonicity.

`DiffusionResult` reports actual accepted midpoint substeps, rejected trials,
maximum accepted error estimate and substep budget/linear residuals. Trial
budget or timestep underflow raises without committing a partial plasma step;
completion on the final permitted trial succeeds normally.
Prescribed bias interpolates linearly between the supplied endpoints; nonlinear
source resolution remains an outer-time responsibility. The moving Strang
solver calls this component at fixed time/bias in each resistive half-step.
