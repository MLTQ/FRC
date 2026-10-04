# Prescribed-current species evolution

`FrozenSpeciesTransport` advances one actual point specific entropy with RK4,
holding the existing point density and Ampere current fixed. The spatial rate
uses the same RegularSpecies skew operator; the physical log-entropy and
quadratic entropy integral are spatial invariants, while RK4 time error is
measured separately. This does not conserve nonlinear internal energy, because
its matching pressure-induced current/magnetic evolution is deliberately absent.

The finite entropy and positive representable EOS energy are checked at every
stage and endpoint. An exception leaves the supplied state unchanged. No floor,
limiter, diffusion or global mean/energy repair is present. Species 0/1 select
the exact electron/ion mass fraction, not an artificial common drift speed.

Dependents: independent shear-translation and discontinuous-trace controls.
This is not a complete coupled plasma, equilibrium or formation experiment.
