# Physical kinetic metric and self-consistent fluid checks

Checks independent continuum Euler/adiabatic identities for finite-amplitude
simple waves, translating two-temperature contacts and rotating equilibrium.
Regular momentum reconstruction/curl and kinetic energy are compared with
analytic smooth fields, including the first annulus. Local conservative
momentum/energy accounting and shared continuity are checked separately from
physical-component energy error. Actual short unforced evolution tests refine
against the independent characteristic solution and preserve atomic rejection.
Conservative balance is normalized by the original uncancelled face throughput:
a stationary uniform flux can produce a roundoff-sized derivative, whose own
tiny norm is not a valid backward-error denominator. The nonstationary contact
and simple wave also test nonzero, conservative transport.
