# iec_particles.py

## Purpose
Actual-mass, singly charged kinetic electrons and hydrogenic ions for the
electrostatic PIC model. Neither electron density nor electron temperature is
prescribed after loading; the Poisson solve determines the charge-separation field.

## Components and contracts
- `Particles`: positions, velocities, tetrahedra, species (0=e,1=ion), common
  physical-particles-per-macro weight and actual ion mass. Charge and mass use
  the same weight; kinetic energy is physical joules.
- `load_plasma`: independent random 3D species positions within a specified
  spherical seed, Maxwellian component temperatures. Globally equal initial
  counts are a declared initial condition, not enforced later neutrality.
- `WallEmitter`: specified uniform-area electron source on grounded wall faces,
  monoenergetic inward surface normals. Poisson emission count has the specified
  mean current. Actual injected charge/energy are measured, not mean-substituted.
  The 1e-9 fraction displacement toward the neighbouring cell is roundoff scale.

## Limits
The source is an idealized emitter, not a model of the user's air discharge.
Wire/wall absorption is handled by the field mesh; no secondary emission,
neutral/coulomb collisions, ionization, fusion reactions, or mass scaling.
Finite macroparticle shape/noise and initial sampling require convergence.
