# Original species physical flux and energy-rate audit

Reconstructs actual n,v,s, h, J,E,B using independent samplers. Electron/ion
velocities are v-beta*J/(e*n), v+alpha*J/(e*n), with physical masses. Total
energy flux is the sum of each species velocity times its kinetic energy plus
enthalpy, plus full E cross B/mu0. It is evaluated directly from the parent
species expression, without a fitted discrete flux or a global correction.

Returns channels: particle number, current, electron log entropy, ion log
entropy, total energy. Current has no evolved charge variable in this
quasineutral model; its reported volume time derivative is explicitly zero,
and actual strong divergence/normal jumps are separate diagnostics.

Energy rates use actual conservative n_t,C_t, physical v_t,J_t,B_t and the
original thermal derivatives. Harmonic constant bias remains in B dot B_t and
Poynting flux. Kinetic, separate internal, current and magnetic powers remain
separate. No diffusion, density floor, heat or energy repair is introduced.
These are independent diagnostics of the full equations, not a solver change.

Restoring the harmonic bias uses a private sampled copy, preserving inputs
even when an independent test sampler returns an array view.
