# Independent physical scalar/field sampler

Builds root-product nodal/edge field bases and cumulative-primitive density
bases, independently of production Vandermonde/moment quadrature. Samples saved
cochains on an independently chosen physical cylindrical Gauss grid. Scalar
coefficients remain nodal; density coefficients remain true conservative means.

Exposes continuous physical L2, RMS, sampled maximum and first-annulus norms,
plus independently integrated cell averages. A zero continuum target reports
null relative error with absolute RMS/max instead of inventing a denominator.
No production geometry, density map, pressure or energy routines are called.
