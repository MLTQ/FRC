# Moving magnetized-fluid trajectory

Actual nonlinear ideal-MHD motion of a compact magnetized sphere in exact vacuum,
using the material fluid and full-vector magnetic force/torque. R=35mm,
Ar+ mass, nc=1e24/m3, central Te=Ti=1eV, Bp=0.3T, Bt=0.15T;
duration R/sqrt(5e(Te+Ti)/(mi+me)). These are a mathematical ideal-MHD fixture,
not a proposed fill pressure, experiment or validated dilute-edge closure.

The initial spherical continuum reference is independently differentiated:
psi=Bp*r²*h^4/2, Btheta=Bt*r*h³/R, h=1-(r²+z²)/R².
Analytic B, curl(B) and JxB plus pressure give initial field/current/acceleration
comparisons. A separate high-order spherical integral gives reference magnetic
energy. The production field uses the actual polygonal material labels, so
initial geometry errors are measured rather than normalized away.

During complete trajectories the all-space ideal-MHD scalar virial relation is
I''=4K+6(gamma-1)U+2Wb, I=integral rho |x|². For gamma=5/3 and conserved E,
I(t)=I0+I0dot*t+2E0*t²-2 integral_0^t (t-s)Wb(s)ds. An accumulated trapezoidal
integral is an independent time-sampled observer, never a force or repair.
The diagnostic includes rotational kinetic energy. Spatial refinement compares
actual boundary/field/current/velocity and magnetic/thermal/kinetic histories.
Energy and virial identities alone do not establish field accuracy.

Current reconstruction gives the within-element curl; Q1 geometry can also
introduce tangential B jumps and distributional face currents. Do not mistake
this sampled current for a validated finite-mass current space. The model still
omits current inertia/Hall/relative species transport/resistance and coil fields.
Reported ion skin-depth and drift scales identify this ideal limit's limitations;
the mathematically cold dilute edge is not assumed physically valid.

Every accepted state and frozen source are preserved. The existing resource
and immutable-source helpers are reused, with no new verification framework,
unit tests, artificial atmosphere, numerical diffusion or energy corrections.

The first r4 trajectory exposed the strong-curl limitation. Subsequent trajectories
retain it separately and report a continuous weak-Ampere projection including
face contributions, with independent initial-current error under refinement.
Frozen initial run is preserved. Dynamics and magnetic forces are unchanged.
