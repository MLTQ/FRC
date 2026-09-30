# Independent physical fluid norms

Uses independent root-product scalar and primitive density reconstruction.
Regular radial/angular momenta use a Legendre-moment fit with independent
20-point integration, not production monomial-fit matrices. These reconstruct
the declared discrete fields; accuracy is measured against separate exact
continuum solutions. Physical volume, first annulus and sampled maximum errors
are retained, with density-signal errors distinguished from full-density errors.

`measure` also integrates true kinetic/two-species energy and physical entropy.
`rates` differentiates their physical expressions using actual n_t, m_t, s_t.
It never imports production thermal or kinetic energy routines. The caller
compares the resulting local integrals at higher quadrature with production.
