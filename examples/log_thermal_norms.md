# Independent physical log-thermal norms

Samples saved log entropy, its rate, conservative n/n_t and E/J/B coefficients
with root-product fields and primitive density at independent Gauss orders.
Computes physical pressure/temperature, specific-entropy rate, entropy-proxy
rate and separate local entropy/density contributions to internal-energy change.
The old lumped proxy rate is sampled with the same nodal basis for a matched
continuous physical norm, rather than comparing different error definitions.

Both species retain physical entropy and proxy integral defects separately.
Nodal specific-entropy rates use the same analytic target for old and new
methods and expose whole/first-four/interior/wall regions. The first-four nodal
region shrinks on refinement; the continuous first-annulus norm is distinct.
Relative error is null for exact-zero targets; absolute RMS/max remain visible.

Independent joint magnetic/current energy rates retain local bias cross work
and the true supplied density derivative. The absolute joint B_t norm is a
spurious-battery measure only in declared curl-free barotropic/isothermal cases;
it is not an error norm for general nonbarotropic pressure evolution.
