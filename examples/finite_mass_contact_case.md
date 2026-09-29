# Contact projection and measured errors

Initial density/momentum use independent exact cell averages. Electron pressure
is the analytic node value; cell ion pressure complements its actual Q1 average
to keep the discrete total pressure constant. This differs from exact continuum
pi/energy means by measured projection error, which must be reported at time
zero as well as the final time. No evolved state is reset to that projection.

The full energy convention is initialized once by the existing state constructor.
Exact continuous contact density, momentum, energy and pi cell means use separate
Gauss quadrature. Electron errors use their nodal representation; they must not
be labeled continuous field norms. Total-pressure/velocity errors expose
spurious dynamics caused by inconsistent nonlinear species transport.
