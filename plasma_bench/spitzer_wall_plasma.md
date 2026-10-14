# Temperature-dependent electrical resistance

`SpitzerWallPlasma` adds conditional Z=1 weak-field scalar resistance to the
complete moving model, retaining electron-ion thermal exchange. Coefficients
are recalculated from each actual nonlinear stage; no prescribed temperature
trajectory enters the simulation. Both species require gamma=5/3.

For physical quadrature n and Te in eV, C=e^4/(4 pi epsilon0)^2:

    nu_ei = 4 sqrt(2 pi) n C lnLambda / [3 sqrt(me) (e Te)^(3/2)]
    eta = 0.51 me nu_ei/(n e^2), Q_J = eta |J|^2.

This uses the rounded 0.51 friction factor and electron heating convention
in the [NRL Plasma Formulary 2023, pp.36–38](https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf).
Its reciprocal differs by 0.04% from the older diagnostic's rounded 1.96.
Joule heat goes to electrons at leading order in me/mi; thermal equilibration
transfers heat between species. This is not a complete finite-mass collision
operator or Braginskii model. Thermal forces, frictional heat flux, conduction,
viscosity, neutrals, ionization and radiation are absent.

`spitzer_samples` reconstructs n, C_s/n and the EOS on electromagnetic
quadrature, using the actual EOS reference density. Existing Coulomb-log
branches validate positive states and slow thermal ions, with no floor.
`SpitzerTransport` retains the full current transport, adding P.T W eta J
to Ohm forcing. Its electron source is (S.T W Q_J)/(V h_e). The two use the
same eta and current, pairing native field/current loss with thermal gain.
The high-order adjoint can yield negative nodal sources even for positive
physical heat; no clipping or global energy correction is applied.

The parent is initialized with zero prescribed resistance and its transport
operator is replaced, not superposed. Parent exchange remains active and
native power diagnostics include the actual variable-resistance heat. The
parent scalar-resistance diagnostic is removed and replaced by actual ranges.
Geometry caches rebind the current density model on every evaluation.

Diagnostics expose actual Hall parameter, Coulomb log and drift relative to
each species' thermal speed. Values above0.1 for Hall/drift request review;
their absence does not prove locality, Maxwellian distributions or validity.
This scalar closure must not be applied as magnetized transverse resistivity.
The mathematical low-field Bessel example measures local scale ordering.

Periodic axial, impermeable radial, fixed psi and zero wall Btheta boundaries
are unchanged. There is no vacuum or formation model. Original solver bounds
remain; no added numerical damping, filter, state floor or retry is included.
See examples/spitzer_bessel_decay.md for the end-to-end physical comparison.

Four complete12us decay runs establish the coupled response and energy-time
convergence; see docs/SPITZER_DECAY_RESULTS.md. Use LiftedModalGauss via
examples/lifted_spitzer_decay for coordinate addition: the original160-step
run rejected floating wall-trace cancellation and remains archived. Current
error does not mesh/time-converge here, despite small absolute errors; do not
infer accuracy of the separate long-contact problem or formation.
