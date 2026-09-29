# Conservative electron enthalpy on primal fluid faces

`scalar_current_flux` reconstructs an RT0 (normal-continuous) integrated current
inside each dual rectangle in x=r²/2,z and integrates its product with a nodal
Q1 scalar through primal faces. It is not a sampled weak Hall current treated
as if pointwise divergence free. Each primal face crosses two half dual cells;
the scalar's exact half-segment means are 3/4 and 1/4 of endpoint values.
Nonuniform radial dual widths, half axis/wall support and periodic z are explicit.
With scalar=1, this maps a divergence-free dual current to a divergence-free
primal current, including the physical radial boundary trace.

`enthalpy_flux` sets scalar=pe/n [J] and multiplies by -gamma/((gamma-1)e),
giving integrated outward watts. This is the physical -(epsilon_e+pe) J/(ne)
term missing from a common-velocity total energy flux. It is separate from
pressure Poynting transport; neither term may be silently dropped because the
pressure electric field has zero curl. No local energy reset is performed.

This interpolation is a candidate second-order flux, not an exact local work
identity. Its conservative divergence, physical face integrals and local
ion-energy defect must be measured independently before full warm coupling.
