# Independent species-rate fixture

Uses the saved analytic magnetic/current/density fixture but computes pressure
by differentiating p=e*n*exp(s)*(n/n0)^(gamma-1) directly. The original current
entropy rate is -a*J dot grad(s)/e. Electron/ion entropy profiles differ and
both density and temperature vary radially and axially. No production gradient,
transpose, pressure-work identity or repaired target enters the oracle.

`initialize` sets the actual point entropy densities. `exact` accepts physical
points, EOS parameters and physical mass fractions; it returns both species'
vector electric pressure fields and scalar entropy-density rates. These are
instantaneous manufactured verification fields, not an equilibrium or formation
experiment. The analytic oracle applies to smooth traces; discontinuous-trace
conservation and transport require separate controls.
