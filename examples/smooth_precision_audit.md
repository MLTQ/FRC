# Original contact induction rejection and corrected representation

Reproduces the exact original 8×16 physical contact and source at the unchanged
linear tolerance, retaining the rejection text, original source fingerprint,
full state/forcing and new split fields. The older production code is unchanged.

Independently replays original Ohm/Ampere and difference-first Faraday equations.
Their backward errors and the full-field rate-relative error are both recorded.
The separable kernel is exactly annihilated before metric division. A large
full-rate replay error remains visible when recombining large/small fields in
double precision; no claim is made that this recombination is well conditioned.

This is an accounting and failure-preservation audit, not an independent exact
solution for the roundoff-level contact magnetic rate. Weak-signal mixed-system
tests and continuum wave/temperature evolution provide separate accuracy checks.
