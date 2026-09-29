# Compatible local nodal and edge polynomials

`bases` returns sparse nodal interpolation N, its derivative, interval-average
histopolation E, and the metric incidence D=(f[i+1]−f[i])/h. On each interval,
an odd-degree local Lagrange polynomial includes both endpoints. Its edge basis
is minus the cumulative nodal derivatives times the corresponding interval
length. Thus dN=E D and each edge basis integrates to its designated interval's
length, with zero integral on other intervals. No physical solution is used.

Radial nonuniform x=r²/2 nodes use physical interior stencils; uniform periodic
z wraps unrolled stencils without a duplicate endpoint coefficient. Degrees
1/3/5 are explicit. This is a local overlapping polynomial construction, not
a claim of a standard spectral-element stability theorem. N is continuous,
its derivative/edge fields may jump at cell boundaries. No slope limiter or
positive-interpolation guarantee exists. Independent polynomial, interval
integral, commuting derivative and physical field tests must verify its use.
