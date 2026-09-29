# Moving-mesh magnetic vacuum and force derivative

`MovingField.solve` prescribes induced a=psi/r² only at material nodes and solves
the stationary weighted magnetic equations at all exterior nodes and multipole
amplitudes. There is no exterior mass, pressure, resistivity or current variable.
`MovingMagneticState` retains the symmetric operator, full solution, material
loads, analytic geometric energy gradient, vacuum residual and two independent
magnetic energy quadratures.

For regular induced fields decaying at infinity,
`Eind=(pi/mu0) integral r³ |grad a|² dr dz` over the full meridional half-plane.
This equals physical `B²/(2mu0)` energy because the cross term is a boundary
divergence whose internal traces cancel. `fields` computes Br=-r a_z and
Bz=2a+r a_r, with optional uniform external bias. Physical finite-domain B²
quadrature plus exact exterior-mode energy independently checks that identity.

The geometry derivative uses the envelope theorem: stationary vacuum values
need not be differentiated, but all material/exterior element geometry terms
must be included. The declared exterior mesh-motion transpose returns those
terms to material coordinates. The material coefficient's ideal-advection and
external-source dependence is additional and belongs to `moving_advection`
and the time integrator. No energy residual is fed back as a force correction.
