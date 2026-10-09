# Fixed radial-wavelength full-RHS consistency study

Selects one independently derived real continuum mode, default radial16,N128,
and holds it fixed while varying fluid elements/points and field radial/axial
cells. Both signed perturbations use the complete original physical equations
and strict induction acceptance, with only the previously verified uniform
physical current-inertia preconditioner. No physical equation is substituted.

Existing independent quadrature16/20 rate reconstruction checks all density,
velocity, species entropy, B, current rate and full E components, initial
representation and axis/wall/max norms. New independent current/face diagnostics
retain strong versus weak and both one-sided traces. Saves all state/rate
coefficients, per-face integrals, expanded-field fingerprints and source hashes.

Every case has a terminal status. A rejected physical solve retains its input
and already computed rates and is explicitly failed; it cannot count as a
passing convergence point. Completed raw cases remain available during a batch.
Fresh full replay is required. These instantaneous tests do not establish
nonlinear stability, complete trajectories or laboratory formation.

The former N64 reference's growth is nearly converged but its real-Nyquist
magnetic value and zero-Nyquist derivative differ by0.325% in an independently
sampled current component. The verified N128 reference is the default target;
the original archive remains untouched and available as an explicit control.
