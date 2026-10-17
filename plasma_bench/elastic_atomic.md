# Ground-neutral electron collisions in the atomic source

`ElasticAtomicMoments` adds the explicitly supplied Ar ground-neutral/electron
elastic pair to `VectorizedCollisionalAtomicMoments`. It selects `network.ground[0]`,
not every zero-charge level, because ground-state cross sections do not establish
excited-neutral scattering. Every existing occupied charged pair remains active.
The source is argon-specific; other gases and unspecified neutral pairs are not
silently assigned this table. Masses come from the same physical atomic source.

`TabulatedElasticPair` supplies K0,H0 using the declared stationary-target energy,
finite domain and isotropic finite-drift Maxwellian assumptions. The inherited
opposite-row pair assembly adds force and lab-energy exchange to both species.
Neutrals recoil and heat; there is no infinite bath or duplicated Joule source.
Collisions advance at fixed populations between reaction halves, as documented
in atomic_collisions.md. All moments pass through the existing conservative
moving-support return; no velocity, temperature or number is reset.

`maximum_elastic_outside_probability` tracks the largest relative-velocity
probability outside the table at evaluated collision stages. It is a domain
coverage diagnostic, not an error bound on unprovided cross sections or heat.
The caller records input identity; raw private tables are not archived.

The standalone full relaxation evidence validates the supplied elastic moments
against independent velocity integrals. Moving reaction/collision splitting,
current/thermal resolution, distribution validity, missing ion-neutral/neutral-
neutral scattering, charge exchange, excited-target channels and wall physics
still require evidence. The fixed positive Coulomb logarithm remains explicit.
