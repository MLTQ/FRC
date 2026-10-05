# Complete wall-conforming trajectories

Separate candidate driver follows the frozen implicit trajectory protocol and
uses StrongWallReference/Gauss throughout. Analytic initial approximation,
physical masses, complete RHS and unchanged nonlinear/induction bounds are
explicit. All stages are in the declared wall space. No state repair or added
damping occurs. The duration is the exact independent contact transit or wave
period, with optional explicitly reported fraction for short controls.

Stores initial/final/penultimate states, every acceptance residual and limit,
iteration history, conservation changes and independently sampled field norms.
A failed next step saves the last accepted state and raises; it is not complete.
Continuous energy, full E/temperature, boundary reactions and original curl
convergence require the separate strong_wall_audit. Frozen original trajectories
and source closures are unchanged. These are bounded verification fixtures.
