# Full moving plasma with driven-wall current

Assembles all prior bulk, electromagnetic, current-advection/stress and species
terms using J=J0(B)+chi*y, with y an independent evolved current. Fixed scalar
eta or actual weak-field Z1 Spitzer eta, paired electron Joule heat, and Coulomb
electron-ion temperature exchange remain. Both gamma=5/3. The existing wall
normal reaction does zero work; actual n_t drives the full current inertia.

The extra Ohm equation supplies y_t. A density-dependent Schur elimination calls
the unchanged strict induction solver in its original space. The complete
unreduced Ohm equation is then checked at the same1e-10 bound, including all
physical current rates and restored species potentials. No acceptance bound
is weakened. Old solver files and frozen results remain untouched.

The actual external lift is Eext,theta=-r Bbar_dot/2. Its forcing enters every
Ohm test, including chi. Source power uses the full J plus harmonic field energy.
The returned electric coefficients describe E_ind; add analytic Eext to obtain
physical E. Joule heat uses the same full current as inertia, Lorentz force,
transport and Ohm. There is no compensating heat or numerical diffusion.

This enriches a filled, periodic rigid-cylinder model. Exact current/field
consistency and continuum wall power still require physical convergence. No
vacuum, real plasma edge, sheath, coil geometry, real ends, conductivity tensor,
thermal conduction/forces, viscosity, neutral or ionization model is supplied.
Full simulations against the independent cylindrical baseline evaluate this
physical change; no unit-test or new verification framework is introduced.

Measured full24µs160/320-step runs and their remaining spatial/current limits
are recorded in docs/WALL_CURRENT_RESULTS.md. The finer spatial run is still
active; this is not yet an accepted formation or vacuum-boundary model.
