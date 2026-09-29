# Prescribed-current physical entropy transport

`LogEntropyTransport` evolves the physical log-entropy scalar through the new
number-weighted Galerkin metric, with the supplied compatible current and
conservative density explicitly held fixed. This is a **transport verification
fixture**, not a frozen-density substitute for the full moving solver. No
magnetic, current, bulk or density evolution is implied by its RK4 step.

Every RK4 stage and endpoint checks finite positive physical thermal fields.
The step creates new arrays and rejects without changing accepted entropy or
the immutable prescribed current. It returns worst original scalar residual
and total linear iterations. Time step selection is the caller's responsibility;
there is no hidden retry, diffusion, floor or energy/entropy reset.

The continuum test uses a closed-compatible toroidal field whose axial current
changes sign radially, so there is no net axial current and no nonzero wall
Btheta. Exact radially sheared periodic entropy advection then tests the local
axis rate over finite time. Physical entropy is a linear invariant of the
discrete scalar transport; thermal energy and entropy-proxy conservation are
independent convergence diagnostics, not forcibly imposed invariants.
