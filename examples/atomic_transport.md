# Full spatial species histories

Two independent physical limits exercise transport of the actual moment state.
They are component verification fixtures, not formation experiments. No unit
suite or reporting framework is introduced.

The plasma contact uses equal electron/ion density n/Nref=1+0.2 sin(2πz/L),
common laboratory velocity 100 km/s and constant separate pressures corresponding
to T=0.2 eV at reference density. Exact density translates through one full
1 m period in 10 µs; charge/current and pressure gradients are exactly zero.
Physical electron and argon ion masses are retained. Reactions/collisions are
disabled in this declared inert contact limit. Grids are 32/64/128 cells with
the physical electron acoustic CFL; the velocity frame removes common bulk
energy without changing the lab transport speed.

The neutral case translates through a 1 cm periodic component domain in 40 ns
at 250 km/s while one real supplied Ar I level decays to a stable lower level.
Initial pressure is constant and T=0.02585 eV at reference density. Both level
densities translate with their independent exp(-At) weights, retaining motion and
energy as the lower state is populated. Electron number/energy/current remain
exactly zero. Grids are 16/32/64 cells, with actual reaction/transport splitting.
These lengths/velocities are verification inputs, not laboratory recommendations.

The independently known cell average contains sinc(π dz/L). For the photon
ledger at each fixed point, integrate the translating upper population:

    P/eDelta = (1-exp(-At))
      + a Im[exp(ikz) A (1-exp(-(A+i k u)t))/(A+i k u)].

Here a includes the cell-average sinc factor. This is cumulative local photon
production, not propagated radiation or a photon field advected with gas. The
reference does not call the numerical flux/source step to construct a solution.

Each run saves every accepted spatial state, time, actual fixed frame, analytic
final state and integral budgets. Invalid stages are preserved as failures,
without normalization, charge/energy fixes or a matter floor. Source/input
provenance is frozen; raw personal atomic tables/rate matrices are not bundled.
Finite numerical diffusion, arbitrary charged-field feedback, radial geometry,
vacuum-front accuracy and physical ends remain unfinished.
