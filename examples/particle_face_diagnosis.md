# Independent normal-face instability diagnosis

Compares the original arbitrary-point two-stencil trace with a single centered
normal face. Builds the complete periodic linear particle-advection operator,
reports every-mode growth and independent RK4 full-transit density-wave error,
and retains spectra/operators/snapshots. Both flow signs, degree 3/5 and three
resolutions prevent a one-direction dissipative substitute from passing.

No thermodynamic solver participates: this isolates the density instability
from pressure floors, entropy projection, Hall terms or a thermal time-step
failure. Nonpositive diagnostic densities remain recorded instead of repaired.
The added source/companion archive supplements the frozen bulk source identity.
