# Next species-pressure and entropy contract — 2026-10-04

This is a derived contract and measured diagnostic, not an implemented species
integrator. Continue frc-7mm.4. The current-transport component is independently
verified; complete thermodynamic and moving evolution are still open.

## Original physical identities

With total mass M, alpha=me/M, beta=mi/M and w=J/(e n),

    vi = v + alpha w,     ve = v - beta w.

The common quasineutral continuity law n_t+div(nv)=0 also equals each species
continuity only when div(J)=0 with consistent boundary/interface fluxes. For
collisionless adiabatic entropy variables Cs=n ss, the original laws are
Cs_t+div(Cs vs)=0. Pressure forcing in Ohm's law is

    E_pressure = [alpha grad(pi) - beta grad(pe)]/(e n).

This follows by taking alpha times ion momentum minus beta times electron
momentum; it retains the physical electron mass. For species gamma and
h=e exp(s)(n/n0)^(gamma-1)/(gamma-1), p=(gamma-1)n h,

    grad(p)/n = grad(gamma h) - h grad(s).

Let a_i=alpha and a_e=-beta. The drift part of advective entropy evolution
is Cs_t=-a_s/e J·grad(s_s). At fixed shared density its thermal power is
-a_s/e integral h_s J·grad(s_s). Its difference from current-pressure work is

    integral J·E_pressure,s - thermal_power,s
        = a_s/e integral J·grad(gamma_s h_s).

This vanishes for a divergence-free current with zero normal boundary flux,
including shared internal face contributions. Strong div(J) and normal jumps
are not exactly zero in the present physical reconstruction. A discrete curl
current annihilates a *compatible incidence gradient* exactly; an arbitrary
nonlinear interpolated enthalpy gradient is not necessarily in that space.
Likewise, global species-entropy conservation requires the corresponding
entropy-gradient/charge-flux contract. Subtracting a global power remainder
would not establish either local equation.

## Saved diagnostic and measured scope

`examples/species_charge_probe.py` takes the verified degree-five saved currents
at 4/8/16/32 fluid elements and prescribes distinct smooth, nonuniform electron
and ion entropies at the actual regular fluid points. Independent reconstruction
at orders 16 and 20 compares actual nonlinear thermodynamic gradients with
compatible incidence gradients of the analytic nodal potentials. The extra
entropy amplitudes and all EOS constants are explicit in the frozen code.
No state is evolved, and no force/temperature/entropy is repaired.

Frozen probe source:
`c6a23e1abb844963d2eababc68a6b8c61a3d808ec8321dbee609b9bec2351c59`.
Evidence: `outputs/species_charge_probe/` with raw input hashes and source archive.
The auditor separately records its source hash and replays both quadrature
orders before checking the work identity and potential-error convergence.

Electron pressure/thermal mismatch using the actual nonlinear gradient is
1.01843e-6, 1.19174e-9, 1.21808e-12, -2.42861e-15 W as the mesh refines.
The compatible potential-gradient work is roundoff-scale in this fixture.
But the compatible replacement also introduces a measured relative pressure
field change: electron 2.95114e-2, 1.61918e-3, 6.03154e-5, 1.97502e-6;
ion 4.36450e-2, 2.34862e-3, 8.74164e-5, 2.86635e-6. Hence exact global work
is possible while finite local interpolation error remains. Strong entropy
integrals happen to cancel to roundoff in this smooth fixture; this is not a
proof for general states, interfaces or nonlinear evolved fields.

## Concrete next implementation gate

Derive a common particle/current flux and species entropy update on the actual
regular state. It must preserve both species entropy integrals and pair each
pressure-current term with the actual grid thermal gradient h_s, while converging
to the original local pressure and entropy rates. Establish the relation to
compatible scalar potentials explicitly, with zero numerical wall flux and
paired internal/periodic faces. Do not substitute another density continuity,
coarse cell-mean temperatures, global energy repairs or a conducting atmosphere.

First test nonuniform density and independent species-temperature profiles,
including discontinuous element traces, against analytic original parent rates;
measure local/axis/wall force and entropy errors, physical versus grid energy,
charge/particle continuity and higher quadrature. Only then combine with the
verified current transport and run unforced complete-parent wave/contact
trajectories, enthalpy/Poynting boundary balances and convergence controls.
Physical heated moving vacuum and formation remain later gates (frc-7mm.5).

## Candidate discrete construction for the next implementation

The following construction is derived but **not implemented or verified**.
Let S be the existing regular scalar sampling matrix, W physical field weights,
V the fluid nodal volumes, P compatible electric sampling and G the compatible
potential gradient. Obtain a common potential phi_s by sampling the actual
entropy ratio (S Cs)/(S n) at field nodes, with explicit central traces at fluid
interfaces; a constant entropy must remain constant. A second potential
phi_gamma_h samples gamma times the actual EOS enthalpy there. Positive
reconstructed density must be checked, not floored. Let h_i be the true grid
thermal gradient at each fluid node, not a redefined continuous energy gradient.
A candidate instantaneous species pair is

    Cs_dot = -a_s/(e V) S.T W [Jh · P G phi_s],
    F_pressure,s = a_s/e [M_E G phi_gamma_h
                         - P.T W ((S h_s) P G phi_s)].

Since S*1=1 and the Ampere current annihilates M_E*G, the species entropy sum
and constant-entropy state have zero drift, and

    sum V h_s Cs_dot = J.T F_pressure,s.

This uses the real thermodynamic gradient and an incidence identity, without
subtracting a global energy residual. It leaves the actual common n_t unchanged.
But S h_s differs at finite resolution from h(Sn,SC), and projecting entropy to
common potentials can spread interface variations. Those are approximation
choices to test, not silently adopted field definitions. Prove a local shared
flux interpretation (including discontinuous traces), convergence to the
original local rates and unforced stability before using this pair. Exact
algebraic work and total entropy alone are insufficient. This is the concrete
starting point for the next implementation under frc-7mm.4.
