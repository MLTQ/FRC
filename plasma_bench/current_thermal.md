# Finite-mass current and thermal subsystem

Prescribed positive quasineutral n and constrained barycentric velocity v=0.
Physical masses define alpha=me/(mi+me), beta=mi/(mi+me), delta=beta-alpha,
mu=mi me/(mi+me), w=J/(ne), vi=alpha w and ve=-beta w. The constraint reaction
is mass-proportional and does zero summed species work. This is not the moving
solver: drift stress, variable-density inertia and bulk advection must follow.

Both species entropies are evolved with their actual opposite current velocities.
Their Ohm pressure combination is beta Ep_e-alpha Ep_i. The current energy is
integral mu J²/(2n e²); `star_field` includes B*=B+curl(mu J/(n e²)) so Hall
transport retains full convective current inertia, not only a time regularizer.
All products use the compatible electric maps and positive physical quadrature.

The stiff solve removes only explicit curl-free scalar gradients. Pressure
potential is returned for full electric/Poynting reconstruction. `hall` sets
me=0 consistently; it is an explicit control, not an adjustable effective mass.
Density, gammas, wall traces and domain retain ElectronPressure/TraceHallCurl
contracts. Density and its inertia coefficient are immutable after construction
so cached matrices cannot silently become inconsistent. No resistivity,
collisions, floors or total-energy reset is included.
Diagnostics sum physical magnetic/current/species energies independently of a
conservative face-flux account; local mismatch must still be measured.
