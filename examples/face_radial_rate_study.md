# Independent full two-dimensional radial-face rates

Nine degree3/5 mesh and quadrature controls compare the actual state derivative
against the original continuum Euler and two-species pressure equations. The
radial pressure, centrifugal and advection accelerations are measured separately
against independent complex-step analytic derivatives. First-annulus norms and
sampled maxima prevent compensating force errors from hiding in equilibrium.

True physical energy rate, separate conservative account and their local/global
defect are saved. Global energy defect is scaled by initial physical energy
times sound speed/radius; no timestep or state repair makes it vanish. All raw
initial/rate arrays, component face rates and independent physical energy are
hashed and saved incrementally. This instantaneous manufactured test does not
establish long-time stability, shock accuracy or full moving-plasma coupling.
