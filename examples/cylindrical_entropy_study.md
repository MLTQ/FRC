# Cylindrical verification protocol and raw evidence

`rate_case` samples the full 2-D original Euler reference, computes actual
pressure/velocity/energy rates independently, and saves full/axis/max errors.
The radial-order sweep uses fixed 128 axial points to expose the axial error
floor. First-element and pointwise errors are retained beside global norms.
Regular physical integration uses 32/48 radial nodes and 4x/8x axial samples.

`evolve` runs unforced full-transit contacts, one radial sound-crossing of a
rotating equilibrium, or one complete small-amplitude Bessel acoustic period.
The last is a linearized reference, with a tenfold amplitude control. Contact
step halving and reverse propagation are distinct protocol cases. RK4 uses a
declared conservative frequency bound; no limiter, floor, forcing or repair.
Every roughly 10 percent accepted state/history is atomically checkpointed;
the last accepted input and actual failure are preserved on rejection. Terminal
records cannot be overwritten. No automatic restart is implemented here.

Global conserved histories exclude radial momentum, which is not an invariant.
Physical energy drift, regular continuous energy, and its quadrature discrepancy
are distinct. Exact zero momenta use a stated 1e-12 characteristic normalization;
that scaled drift is not a relative error in a nonzero invariant. Specific
entropy change is not an error for advected nonuniform contact entropy; for
uniform-entropy acoustic modes it measures spurious local entropy evolution.

The source closure is frozen before any data are written. Full arrays, physical
volumes and time snapshots support independent audit and original-platform
replay. These are positive-density fluid fixtures, not formation experiments.
