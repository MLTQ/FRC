# Curl-kernel precision separation

The unforced translating contact exposes an ill-conditioned electric background:
Ez(z) is independent of r and produces no curl, while roundoff-sized remainder
coefficients carry the magnetic rate. Reconstructing their sum in double precision
cannot guarantee a 1e-10 error relative to that tiny magnetic rate. The original
solver rejects despite original full-equation backward errors near roundoff.

`separate` algebraically extracts Er/r independent of z and Ez independent of r,
including a possible periodic harmonic axial voltage. These are exact members
of the declared Faraday kernel, not a small-field cutoff or fitted correction.
No physical field is discarded. The complete electric field is restored for
Poynting transport, and the unchanged inertia/current equations act on the curl.
`difference_curl` subtracts equal coefficients before metric division so those
kernel fields vanish exactly. It is the same mathematical discrete curl.

`KernelInduction` retains all original strict relative residual checks for the
non-kernel solve at the caller's unchanged tolerance. Original full Ohm/Faraday
and reconstruction errors are additionally checked as backward errors using
the original terms' magnitudes. Full-E Faraday error relative to the tiny rate
is reported separately and can be large when recombination is ill-conditioned.
This is an explicit numerical representation change, not a claim that a tiny
rate can be recovered accurately from a rounded large full electric array.

Mass inversion and subtraction can leave roundoff-level forcing noise. Independent
mixed-system, kernel addition/weak-signal and contact tests must bound its effect.
The subset of separable kernel fields is not a general Helmholtz decomposition;
an arbitrary nearly curl-free field may still fail. Failures remain atomic.
The historical InertialInduction/MovingCurrent and their sources are unchanged.
The smooth reference alone selects this solver; residual diagnostics distinguish
the checked small-signal maximum from descriptive conditioning ratios.
An exactly zero remainder has a null conditioning ratio rather than an invented
finite denominator or a non-serializable infinity; zero forcing remains zero.
