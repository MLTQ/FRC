# Driving and momentum refinement across a rejected solve

Independently reconstructs original bulk/Lorentz/transport/species/reaction
momentum terms and wall-projected negative curls of the ideal/Hall/pressure/
transport driving electric fields. Compares fixed-reference original parent
terms at quadrature16/20, retaining every physical radial-cell squared norm
and volume/axis/wall errors. Uses bounded radial chunks without changing fields.

The driving sum excludes inertia because a rejected induction solve supplies
no accepted Jt. Its analytic target is the corresponding ideal+Hall+pressure
sum. Error normalization by the full continuum B rate is explicitly separate
from own-driving-term normalization. Neither value is an achieved B-rate error.
No failed full case becomes passing; every output retains original status and
full_plasma_rate_accepted=false. No zero or surrogate inertia is introduced.

Check shared accepted-baseline coefficients against the earlier full capture,
and repeat the complete diagnostic for bitwise arrays and metrics. These
forcing diagnostics guide further investigation but do not satisfy complete
induction, trajectory, local-interface or formation acceptance.
