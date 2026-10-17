# Tabulated elastic momentum and heat exchange

`TabulatedElasticPair` supplies finite-drift binary elastic moments for two
isotropic Maxwellians and an explicit ELASTIC momentum-transfer table. It uses
physical masses and retains target recoil. Table energy is stationary-target
projectile energy E=ma*g²/2, not reduced-mass energy. A stationary-heavy-target
electron table is used as a COM momentum-transfer table; that angular convention
is a heavy-target approximation, not a finite-mass scattering recalculation.
EFFECTIVE tables also contain inelastic contributions and are rejected.

The force follows the binary transfer integral in
[Pinto & Galli, section2, equations1–12](https://arxiv.org/abs/0804.2856).
The accompanying energy law below is derived from the same elastic events:
delta(pa)=mu*(g'-g), delta(Ea)=V.delta(pa), |g'|=|g|.
With T in joules, u=Ua-Ub, M=ma+mb, mu=ma*mb/M, S=Ta/ma+Tb/mb,
relative velocity g is Gaussian with mean u and scalar variance S. The
conditional pair-centre mean is Vcm+(Ta-Tb)*(g-u)/(M*S). Thus

    G = <g_vector*|g|*sigma_mt>, B = <|g|³*sigma_mt> - u.G
    K0 = mu*(G.u)/|u|², H0 = mu*B/(M*S)
    Ra = -na*nb*K0*u = -Rb
    dEa = Ra.Vcm + na*nb*H0*(Tb-Ta) = -dEb
    Qa = (mb/M)*na*nb*K0*|u|² + na*nb*H0*(Tb-Ta)
    Qb = (ma/M)*na*nb*K0*|u|² - na*nb*H0*(Tb-Ta).

No separate prescribed drag or thermal exchange should be added to these terms.
Signed heat is retained. This moment operator leaves populations/internal levels
unchanged. Ionization and excitation require separate reaction channels. For a
future symmetric ion/parent-atom table, its momentum-transfer cross section may
already include resonant charge exchange: separate charge-exchange drag must not
be added to that same contribution. The current electron/Ar table does not
supply ion/neutral exchange.

Angular integration is analytic. For r=g/sqrt(S), d=|u|/sqrt(S), a=r*d,
the radial probability is r²*exp(-(r-d)²/2)*(1-exp(-2a))/a/sqrt(2*pi).
Its conditional mean cosine is L(a)=coth(a)-1/a. A small-a series evaluates
L/a without cancellation; zero drift uses its exact1/3 limit. Radial positive
Gauss–Legendre integration covers every supplied segment, interpolating sigma
linearly in energy. Additional knots around the shifted Gaussian resolve narrow
peaks; no domain is removed at12 thermal speeds. Both cold/no drift gives zero
elastic exchange; a cold moving pair uses the table at its actual relative speed.

Only the supplied energy domain enters. `probability_outside_table` reports the
noncentral-Maxwell probability outside it. This is NOT a bound on missing force
or heat: unknown cross sections and velocity weights matter. No extrapolation,
population/temperature floor, artificial neutral atmosphere or numerical
diffusion is added. Discretization errors are table interpolation and quadrature;
drifting Maxwellian shape, anisotropy and finite-domain truncation remain limits.

`examples.elastic_pair_relaxation` performs full source evolutions and compares
against direct radial/polar velocity integrals and an independent constant-cross-
section zero-drift formula. `ElasticAtomicMoments` supplies its explicit
ground-neutral/electron moving-source coupling. Neither component establishes
a self-consistent distribution or a laboratory prediction.
Private provider tables are explicit local inputs and are never bundled.
