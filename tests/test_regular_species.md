# Species coupling verification

Checks original independently differentiated pressure/entropy fixture convergence,
both physical species' entropy, frozen quadratic measure and actual grid thermal
work. Random discontinuous element traces exercise the common projection and
transpose beyond the smooth fixture. A complete small frozen operator verifies
the number-weighted skew metric and its imaginary spectrum. Invalid density,
EOS overflow and invalid steps must fail with unchanged input.
The three-element discontinuous fixture explicitly uses twelve radial field
cells so field cells subdivide fluid elements as required by the physical map.
Uniform but unequal species temperatures at varying density exercise the exact
pressure curl-kernel separation, independently of work conservation.
An integer-state regression checks that returned physical rates are evaluated
in float64 rather than silently truncated by the input container dtype.

Passing invariants alone is insufficient: the point-interpolation skew control
is deliberately required to retain its large original local-rate error. Full
continuous norms, temporal evolution and physical/mesh/quadrature controls are
recorded by the separate study protocol; no formation claim is made here.
