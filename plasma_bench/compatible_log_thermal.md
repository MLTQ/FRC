# Physical log-entropy pressure/work component

The nodal state is s=ln(w/e), w=p/[n(n/nref)^(gamma−1)], with physical entropy
proportional to integral n s dV. Sampling s with the compatible scalar basis N
and **the same reconstructed conservative n as current energy** gives
p=n e exp(Ns+(gamma−1)ln(n/nref)), epsilon=p/(gamma−1). Exponentials preserve
positivity; overflow, underflow to zero and nonfinite states reject, never clamp.

For electron-oriented current transport, n s_t=J dot grad s/e. Galerkin projection
uses M_n=N.T W n N, so physical logarithmic entropy rate integrates to
integral J dot grad(Ns)/e=0 for Ampere current and the closed compatible boundary.
No lumped graph rate is used. The proxy integral of n w is a separate diagnostic:
it need not be exactly conserved by this finite-dimensional log-state transport.
Its continuum error must converge; no shock entropy inequality is claimed.

Let h=M_n^-1 N.T W epsilon be the number-weighted projection of epsilon/n.
The pressure field is E_s=(h grad s−gamma grad h)/e, with the product projected
onto the compatible electric basis. Continuously this equals −grad p/(ne).
The same projected h gives dH/dt=J.T M E_s. A constant h anchor separates
E_reduced=(h−h0)grad s/e and phi=(h0 s−gamma h)/e, preserving full E while
removing an exact gradient from the stiff solve. No independent pressure fit is
used. Physical electrons multiply this pair by beta and ions by −alpha.

For actual caller-supplied n_t, the full internal-energy derivative is
epsilon_t=epsilon (N s_t+gamma n_t/n). Its density pullback is
R.T W gamma epsilon/n. Local pressure, energy, physical entropy and proxy rates
remain explicit. The current-driven pair supplies only one part of s_t; bulk
advection and its shared continuity/pressure work still must be coupled. These
derivatives do not supply missing bulk energy or reset a conserved state.

The scalar solve's residuals are returned. Current, thermal and geometry use one
density-map object. Boundaries remain periodic z, regular axis/fixed-flux wall and
natural Btheta(R)=0. No vacuum, collisions, diffusion, limiter, time integrator
or complete two-fluid conservation theorem is supplied by this component.
The returned `energy_per_particle` is h, not enthalpy gamma h. Current coefficient
shape/finiteness is checked, and the cached density power is immutable.
