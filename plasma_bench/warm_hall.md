# Coupled warm Hall reference

`WarmHallMHD` composes ideal fluid/species half step, pressure/current half step,
Hall full step, pressure/current half step and ideal half step. Hall density is
the positive Q1 sample of actual overlapping ion mass at the nonideal stage.
Pressure and Hall substeps hold n and momentum fixed. Pressure stages use the
existing magnetic-midpoint/logarithmic thermal mean, with their full electric
field in Poynting work. Total energy changes **only** through face fluxes;
it is never reset to the computed electron/magnetic change. Local numerical
ion-energy residuals remain and must converge.

The ideal half steps use shared conservative bulk mass/entropy transport.
Their existing entropy/ion partition limitations remain. The palindromic
composition has second-order consistency in smooth cases but ideal RK2 is not
exactly time reversible; actual wave/time/boundary controls
must establish accuracy. The fluid CFL is not a sufficient whistler phase or
pressure Picard convergence criterion; `advance(dt_max=...)` allows explicit
accuracy caps, and failures raise without hidden retries or floors.

Every one of the five substeps must validate before accepted fields, U, time,
counts or diagnostic maxima change. The Hall computational preconditioner may
change on rejection, as in the existing integrator. Both nonlinear iteration
histories and residual maxima are exposed; no physical Hall heat is added.

Bounded cylinder: periodic z, impermeable slip fluid, fixed psi at axis/wall,
zero physical Btheta(R) and free Ez(R). This is the explicit magnetic trace
boundary, not a conducting electric wall or a vacuum interface. Different
species gammas, finite electron mass, resistivity/heating partition, collisions,
conduction/ionization and formation coupling are unsupported. Existing cold
Hall and ideal species paths remain available as separate limits/baselines.

`species_rhs(state,U)` returns the full coupled instantaneous fluid derivative,
Udot, zero heat and ideal reconstruction statistics. `rhs` also requires the
matching explicit entropy argument and returns derivative/heat; it raises if
U is omitted instead of silently returning only ideal or mismatched dynamics.
`step`/`advance` use the palindromic composition with implicit nonideal stages.

The fine pressure-only control now rejects a negative ion partition, and its
continuum frozen symbol has a nonhyperbolic direction on the original reference.
Thus the composition is experimental: no uniform stability claim follows from
implicit Hall integration. See WARM_PRESSURE_SPLIT_LIMITATION before extending
it; joint Hall/pressure integration must be investigated without energy repairs.
