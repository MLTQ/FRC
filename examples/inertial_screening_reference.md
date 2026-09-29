# Collisionless carrier-limited spherical screening

For eta=0, zero initial current/field and frozen uniform density, canonical
flux a+mu0 de² Jtheta/r stays zero. The interior therefore satisfies
`a-de²(a_ss+4 a_s/s)=0`. Continuous dipole matching to uniform bias B imposes
`R a_s(R)+3a(R)=3B/2`. Thus
`a=(3B/2) i1(s/de)/[(s/de) i0(R/de)]`, with modified spherical Bessel functions.
Exterior total a is B/2 plus the matched decaying dipole; there is no current
sheet. Current is `-r a/(mu0 de²)` and central Bz is `B/i0(R/de)`.

Scaled modified Bessel functions and a central series avoid exponential overflow
and cancellation. The reference uses no finite-element matrices or solver state.
For de much greater than R, induced current tends to `-ne e² r B/(2 me)`:
it vanishes with carrier density while total field approaches vacuum bias.
The model remains magnetoquasistatic; this is no claim about instantaneous
electromagnetic propagation, large-drift kinetic response or a moving edge.
