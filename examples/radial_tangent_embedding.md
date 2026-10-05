# Exact subdivision of a saved radial tangent

`embed` evaluates each original fluid element polynomial at the nodes of its
integer subdivisions, retaining distinct left/right interface traces. The
transverse first-element polynomial uses the declared degree-four dependent
axis representation; inactive axis slots remain zero. All other components use
the degree-five conservative polynomial. No averaging, filtering or correction.

Requires identical degree/radius/axial grid. In the stationary contact fixture,
base density is radially constant, so this also preserves velocity quotients.
That condition must be checked by callers; this routine is not a generic
nonlinear state transfer. Used only to compare spatial tangent rates on fluid
refinements while retaining the exact original magnetic coefficients.
