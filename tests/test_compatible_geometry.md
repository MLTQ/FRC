# Higher-order compatible physical verification

Checks independent polynomial nodal values/derivatives and interval integrals
on nonuniform and periodic grids. A cylindrical polynomial psi/g fixture has
known B and curl(B)/mu0, respects the actual mixed wall trace, and independently
checks physical energy and magnetic/electric work. A separately specified
electric polynomial verifies original Faraday signs/metrics; a general periodic
potential checks zero curl. Variable-density conservative averages recover a
known physical velocity without fitting it to current. An independent original
three-block system checks the unchanged finite-inertia elimination.

The magnetic/electric work check uses the noncancelled term magnitudes for a
floating-point backward-error bound, rather than dividing by the small net
random work. Analytic local current and field accuracy are checked separately;
this algebraic bound is not a relaxed physical accuracy threshold.

These component tests do not establish nonlinear moving energy/pressure
coupling, shock/entropy robustness or formation. Independent wave rates and
physical quadrature controls provide separate measured convergence evidence.
