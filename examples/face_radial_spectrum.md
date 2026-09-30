# Cylindrical radial sound verification

Linearizes around uniform gas with constant entropy and zero flow. The actual
independent radial-face lift gives the same momentum/continuity field and
physical metric, yielding the energy-paired pressure gradient. Centered
differences of the original nonlinear RHS verify the proposed matrix on its
fundamental acoustic mode at three amplitudes; this is not inferred only from
the energy algebra.

Independent continuum frequencies are c times the first three positive zeros
of J1 divided by the wall radius: regular J0 pressure with zero wall-normal
velocity. Mesh/degree and quadrature controls preserve every matrix/eigenvalue.
An independent matrix exponential checks two-period physical energy and three
RK4 time steps separate time error from the spatial spectrum. The timestep is
bounded using the largest represented frequency. No dissipation is introduced.

This verifies the stationary radial acoustic limit, not nonlinear stability or
an entropy theorem. Source identity and this separately hashed diagnostic are
preserved without altering the frozen trajectory/rate-study closure.
