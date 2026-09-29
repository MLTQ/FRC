# Compressible Hall / resistive fluid integration

HallMHD advances the full toroidal conservative state with explicitly compatible
magnetic energy. It rejects legacy-energy initial states, and requires a finite
positive ion mass. Fully ionized singly charged carriers obey ne=rho/mi;
electron inertia and electron pressure are zero, ion pressure is adiabatic.
The bounded conducting slip cylinder is filled with fluid and periodic in z.
This is not the physical free-vacuum formation or separate-temperature model.

A Strang step is ideal fluid SSPRK2 for dt/2, nonlinear Hall/resistive implicit
midpoint for dt, then ideal fluid SSPRK2 for dt/2. The existing ideal flux and
Maxwell stress supply momentum once. Hall adds no separate fluid force or heat.
Density is frozen only within the magnetic split substep; positive limited
cell reconstruction supplies its actual common-quadrature values. Resistance
and Hall share current, midpoint curl, Poynting flux and positive physical heat.
Total energy receives flux divergence only; heat is diagnosed, not added twice.
The prescribed-resistance closure assigns that heat to the represented gas
internal energy. It does not solve the physical electron/ion heating partition
or justify keeping real resistively heated electrons cold. Those thermal and
transport closures remain required for a physically predictive formation run.

Every trial must pass density/pressure and nonlinear/linear residual checks
before state, time, heat and limiter/iteration totals commit. Input/accepted
states remain intact on failure; the linear preconditioner cache may change.
`rhs` provides the unsplit instantaneous full equations for verification.
`hall=False` disables Hall explicitly while retaining compatible resistivity.

`timestep_limit` bounds the explicit ideal-fluid steps. Implicit midpoint
removes an explicit whistler stability restriction but NOT Hall phase error.
Callers must use `advance(..., dt_max=...)` and measured time convergence to
establish accuracy. Nonlinear convergence or post-Hall positivity can impose
smaller steps. No automatic floor, fitted axis correction or energy reset exists.
Diagnostics count four accepted ideal reconstruction stages per complete step.
