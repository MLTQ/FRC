# Physical log entropy, pressure work and the common density metric

The previous graph pressure bridge conserved a lumped thermal energy but had a
nonconverging local entropy-rate error near the axis. The new component uses the
same physical quadrature, continuous nodal scalar space and **actual conservative
density reconstruction** as the compatible fields. It is a new thermal contract,
not a rescaling of the old graph rate or a correction fitted to the analytic test.
The former implementation and its measured defect remain available.

## Thermodynamic state and metric

Let N sample nodal scalars, R reconstruct true annular cell means, W denote
physical cylindrical quadrature, n=R n_cell and n_t=R n_cell,t. Define

    s = ln(w/e),  w = p/[n (n/nref)^(gamma−1)],
    p = n e exp(Ns) (n/nref)^(gamma−1),
    epsilon = p/(gamma−1),
    S = integral n Ns dV.

S is proportional to physical ideal-gas entropy, up to the constant reference
gauge. Positive pressure follows from the exponential for every representable
finite state. Overflow, underflow to zero and nonfinite values reject; there is
no pressure floor or energy reset. The density is the immutable current model's
map, with no alternate nodal density or independent continuity equation.

The number-weighted scalar mass is M_n=N.T W n N. It is positive definite for
positive quadrature density. A tensor action applies this metric; the unweighted
physical scalar Gram matrix is only a preconditioner. The full original scalar
load residual is checked at the stated tolerance, with iteration counts saved.

## Current-driven entropy and pressure

For the electron-oriented current component, n s_t=J dot grad s/e. Its Galerkin
rate is

    M_n s_t = N.T W [(J/e) dot grad(Ns)].

The gradient is sampled through the compatible incidence. Therefore constant
entropy has exactly zero gradient and zero current-driven rate, including at the
axis. For Ampere current, integration gives

    S_t = integral (J/e) dot grad(Ns) dV = 0,

using M J=−C.T H deltaB and C grad=0 at the closed bounded trace. This is physical
logarithmic entropy conservation, subject to the scalar/current solve's numerical
residuals. The polytropic proxy integral integral n w dV is **not** another exact
invariant of this finite-dimensional log-state approximation. Its signed defect
and convergence are measured separately. No shock entropy inequality is claimed.

Thermal energy H_s=integral epsilon dV has derivative H_s,s=N.T W epsilon.
Define its conjugate h by

    M_n h = N.T W epsilon.

Thus h is the number-weighted projection of internal energy per particle,
epsilon/n; it is not enthalpy (which is gamma h in the continuum). The paired
pressure electric field is

    E_s = [h grad(Ns)−gamma grad(Nh)]/e.

The product is physically projected onto the compatible electric basis. Since
the full gradient does no net current work, the same h supplies

    dH_s/dt = J.T M E_s.

Continuously, epsilon/n=(n/nref)^(gamma−1)e exp(s)/(gamma−1), so the expression
reduces exactly to E_s=−grad p/(ne). Both this continuum identity and the discrete
work pairing must be checked; the latter alone is not an accuracy result.

A constant h0 separates E_reduced=(h−h0)grad s/e from
phi=(h0 s−gamma h)/e. The stiff induction solve receives only the reduced field;
the full electric field restores grad phi. This algebraic gradient separation
retains weak magnetic rates without claiming that the rounded full E array has
the same relative curl accuracy. Original backward and rate-relative errors
remain separate. For finite species masses, electrons multiply their pair by
beta=mi/(mi+me), ions by −alpha=−me/(mi+me).

## Actual changing density and the unfinished moving contract

At supplied n_t, the complete local internal-energy derivative is

    epsilon_t = epsilon [N s_t + gamma n_t/n],
    dH_s/dn_cell = R.T W [gamma epsilon/n].

Physical entropy and proxy rates similarly retain n_t Ns and n_t w. These are
directional derivatives of the actual physical energy/state, independently
verified without repair. The current-driven component supplies only part of
s_t. Full bulk advection must use the same conservative continuity flux, and
its pressure/density work must balance the momentum and total-energy equations.

Combining both species' current pressure rates with the variable-inertia solve
at arbitrary supplied n_t leaves the declared exchange

    d(W_B+K_J+H_e+H_i)/dt
      = integral [mu |J|²/(2ne²)+gamma_e epsilon_e+gamma_i epsilon_i] n_t/n dV.

This must exchange with the missing bulk dynamics; it is not external heating
and is not forced to zero. The earlier cell total-energy/ion-pressure partition
cannot simply be interchanged with two nodal log-entropy fields. The number and
placement of thermal degrees of freedom, local energy flux and conservative
bulk reconstruction need an explicit common contract before full adoption.

## Verification and limitations

Rate tests compare independent physical quadrature, nodal/axis values, both
species' pressure/temperature, log-entropy and proxy rates, local energy and
original induction equations. Constant entropy is a curl-free barotropic
control. Uniform temperature at varying n is generally not exactly representable
by a finite log-entropy polynomial: initial temperature error and false magnetic
battery are measured under refinement, never declared identically zero.

A separate prescribed-current transport fixture evolves a radially sheared
periodic entropy wave over a complete axis transit. Its exact PDE, closed
energy/proxy integrals and independent physical matrix exponential provide
space/time targets. Holding current and density fixed is explicit in that
fixture; it is not a substitute for full moving finite-mass evolution. RK4
validates each stage and returns a new state, without hidden retries or resets.

All boundaries remain periodic z, regular axis/fixed-flux wall and natural
Btheta(R)=0. No numerical atmosphere, vacuum carriers, physical transport,
collisions, ionization, shock stabilization or formation result is supplied.

## Related primary literature

The wider principle of deriving discrete thermodynamic work from a variational
weak formulation is developed for fluids by
[Gawlik and Gay-Balmaz](https://arxiv.org/abs/2211.08745) and for classical MHD by
[Gawlik, Gay-Balmaz and Manach-Pérennou](https://arxiv.org/abs/2503.01239).
Those papers treat broader discrete thermodynamic balance laws. The specific
cylindrical finite-species-mass current/entropy construction above is derived
and verified here; it does not inherit their full-method stability or
fully discrete conservation claims.
