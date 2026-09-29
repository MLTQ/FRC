# Complete moving parent-species manufactured reference

Integrates the analytic reference's original species kinetic/internal energy
and momenta in physical cylindrical volume, independently of production moments
or energy matrices. Analytic derivatives provide exact cell-average conservative
rates. Explicit parent force and heat supply fluid momentum/total-energy,
electron-entropy and electric source terms; no numerical residual is fitted.

Initialization preserves exact averaged density/momenta, exact corner flux,
annular g and prescribed species pressures, but uses the declared discrete
magnetic/current/bulk energy convention once. Its total energy representation
error against continuum remains measurable and is never reset during evolution.
The original nonuniform-density, rotational/compressive, non-Beltrami fixture
and distinct species gammas are retained.

The scalar source at the axis is evaluated at r/R=1e-10 because the independent
cylindrical parent operator contains explicit 1/r. This reference-only analytic
limit must be independently checked against another radius; it is not a carrier
floor, a production field regularization or an exterior atmosphere.
