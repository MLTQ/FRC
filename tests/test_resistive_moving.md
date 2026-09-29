# Moving resistive/source integration checks

Zero eta must reproduce the original ideal trajectory exactly. A nonuniform
material deformation checks inversion of flux labels including the regular
axis. Full evolving FRC/bias data with physical heat verifies second-order
Strang/Verlet time refinement of velocity and response-energy error; cell mass
is unchanged and entropy/physical heat increase without a correction.

An injected failure in the second complete adaptive resistive half-step verifies rollback after
the first half-step and ideal motion have already changed the trial state:
accepted positions, velocities, labels, entropy, heat, time, step count, work
and residual accounting all remain intact. This is failure atomicity, not a
guarantee of mesh validity under arbitrary physical motion.
