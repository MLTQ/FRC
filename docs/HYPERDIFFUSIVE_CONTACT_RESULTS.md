# Moving contact with numerical short-wavelength damping

This is a concrete spatial-solver candidate, not an accepted physical closure.
The collisionless pressure/current growth remains unresolved. The opt-in model
adds sixth-order magnetic/current damping, paired with an explicitly numerical
electron heat source. All original ideal, Hall, inertial and pressure terms
remain; stage/endpoint tolerances are unchanged. Read the companions for
`hyperdiffusive_wall_plasma`, `hyperdiffusive_modal_gauss` and
`hyperdiffusive_contact` for equations, wall conditions and assumptions.

At F8/64, B32/64 the prospective mesh-based strength 0.25 gives
eta6=1.1876780167018035e-19 ohm meter^5. The coefficient is not fitted to the
observed error growth and is distinct from physical resistivity. Its prescribed
electron heat fraction is 1. This is a filled periodic 70 mm diameter, 240 mm
long cylinder with the original 5 mT bias and nonuniform argon contact fixture,
not the proposed laboratory chamber, an FRC or a formation experiment.

## Independent continuum comparison

`outputs/hyperdiffusive_contact_reference_s025` contains a direct matrix
exponential of the original Fourier/Bessel continuum equations with and without
the added analytic eta6 curl^4 J term. It uses the unchanged original growing
mode as its initial condition. At 100 us, the numerical term changes density
by 1.8179997e-7 of its signal, Btheta by 6.6179450e-7 and Bz by 4.8129909e-6.
The original exponential agrees with the independently saved growing mode to
3.91e-11 or better in every component. The calculation exited 0.

These are linear continuum distortion estimates. They do not demonstrate
discrete stability, nonlinear accuracy or small effects on unresolved modes.

## Actual nonlinear evolution — in progress

The unchanged negative 1e-4 perturbation is running for 100 us / 60 steps in
`outputs/hyperdiffusive_contact_f8_b32_s025`, with its existing resource guard
in the adjacent `_control` directory. Actual session: 40062 / PID41757. Latest observed progress:22 accepted steps,
36.666667us, peak ownRSS4.749GB and minimum free disk10.363GB. No automatic
replay or unit tests were added/run. Do not restart the run from this note.

The first five accepted steps cover 8.333333 us. The existing independent
density observer (`outputs/hyperdiffusive_contact_density5`) exited 0; q20
relative density-signal errors at 1.667, 3.333, 5, 6.667 and 8.333 us are
6.99178e-6, 1.32272e-5, 1.38555e-5, 1.25132e-5 and 1.37160e-5.
The next accepted prefix through18.333333us (11 steps) has density error
3.30689999e-5 of the reference signal, versus0.0696199203 for the otherwise
identical original undamped run at the same time (about2105 times lower).
The existing observer in outputs/hyperdiffusive_contact_density11 exited0
(session75938). This is promising intermediate evidence, not a full-duration
pass or a convergence claim. Await actual completion/failure.

The same saved step11 was independently integrated with the existing
matched_mode_physical.changes and matched_mode_energy_scale.measure helpers.
At q20, continuous energy change is -2.0329614e-18 J, or -1.2301322e-7 of the
initial perturbation energy1.6526365e-11 J. q16 gives -1.2301272e-7 on that
same scale. Particle-relative drift is1.0622e-18. The quadratic signal energy
is a normalization, not a conserved invariant of this nonuniform background.
Temperatures and densities remain positive. Exact inputs and helper source
closure are in outputs/hyperdiffusive_contact_energy11; observer24536 exited0,
peak1.571GB. It calls no production RHS or evolution replay.

The next density observer through step19/31.666667us exited0 (session3393).
Its q20 error is6.42403493e-5 of the predicted perturbation. At28.333333us,
error is5.63945143e-5, versus790.2581 for the original undamped state at its
last accepted step. This is a large improvement in the actual trajectory;
full-duration success and vanishing-diffusion convergence remain unestablished.
See outputs/hyperdiffusive_contact_density19 for every accepted time.

## Limits and next action

Gauss is A-stable but not L-stable; extremely stiff decaying modes need not
decay rapidly per time step. The repeated variational curl supplies numerical
higher-order wall conditions. The heat projection is not nodally positive.
No density/pressure floor or stored-state filter is applied. Exact vacuum,
physical coil/end coupling, transport/ionization and formation remain open.

Metadata erratum: the initial source/output scope calls radial cells nonuniform
in r. The actual field grid is uniform in r and nonuniform in x=r^2/2;
h=R/nr is the actual radial field spacing. This wording does not affect the
implemented coefficient or equations. Preserve the running source unchanged.

NEXT collect session 40062, compare the saved trajectory with the existing
continuum observer, and measure actual energy change. A complete accurate
trajectory must precede any claim that this is a useful numerical remedy.
