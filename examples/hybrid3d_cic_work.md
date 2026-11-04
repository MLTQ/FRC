# hybrid3d_cic_work.py

## Purpose
One passive, independent analytic CIC work comparison for a concrete precision
change. Actual full native particles come from the bounded startup replay;
there is no new physical simulation, unit suite or energy correction.

## Components and contracts
`run` reconstructs and checks the recorded device box, then compares CPU64,
original world-coordinate Metal and split-coordinate Metal on exactly the same
positions/velocities/weights. The exact manufactured electrostatic quadrupole
E=(-x,-y,2z), coefficient 1 V/m², has div E=curl E=0. Trilinear shapes reproduce
linear E, so integrated deposited J dot E equals sum(q w v dot E(x)).
Exact total number and first position moment are checked independently too.
Periodic stencil crossing is rejected because this analytic E is not periodic.

Report absolute power error, error relative to signed net power and error
relative to sum of absolute particle powers; cancellation can make net-relative
error misleading. CPU64 is only the regional-array comparison reference, not
the independent analytic work/number reference. Timings include full host
coordinate preparation, conversions, uploads, kernel and collection after JIT.
They are CIC timings, not a complete simulation speedup.

The existing uniform-density/current fixture is also run for both Metal classes,
including ordered, shuffled and partial SIMD-group particles. Source copies,
hashes and native snapshot fingerprint stay with output. No live solver or
running campaign is changed, and this field is never imposed on a physical run.

Usage: `python -m examples.hybrid3d_cic_work NEW_OUT --state NATIVE_NPZ
--setup SETUP_JSON [--dx .01 --dz .01]` with GPU-capable execution.
