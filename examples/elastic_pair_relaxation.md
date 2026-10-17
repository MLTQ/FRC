# Complete elastic-neutral source evolutions

Run `python -m examples.elastic_pair_relaxation FRESH_OUTPUT --cross-sections
PERSONAL_LXCAT_TEXT`. The input must contain exactly one Ar ELASTIC momentum-
transfer table. It is read locally; provider arrays/rate curves are not archived.
SHA, provider attribution and energy bounds accompany the saved physical histories.

Two collision-only problems evolve both velocities and temperatures for eight
initial small-drift relaxation times. A synthetic constant cross section with
equal physical argon masses exposes recoil and equilibration; the sourced BSR
electron/argon table uses electron/39.948u masses, n=1e18/1e20m^-3 and initial
Te/Tn=2/0.02585eV. Both begin at three combined thermal speeds of relative drift
and a nonzero common boost. Populations are fixed; there is no sustaining field,
reaction, spatial transport or chamber model. The equal-mass table spans0–1e5eV;
its constant sigma2e-19m² is a verification fixture, not argon collision data.

Independent references integrate radial and polar binary-event momentum/energy
transfer directly. Opposite polar directions are paired with expm1 to retain
small drift, without the production analytic angular formula. Full8x32 and16x64
radial/polar quadratures are evolved by DOP853 with dimensionless rtol2e-11,
atol2e-13; their complete histories measure reference quadrature uncertainty.
Both integrate the same linear-in-energy table on its entire supplied domain.
The conditional Gaussian centre velocity analytically integrates the remaining
three velocity dimensions. This assumes isotropic Maxwellians, as the production
model does, but uses a separate moment evaluation and time integrator.

Full100/200/400-step RK4 trajectories record every state, reference, momentum
and physical total energy. No thermal correction, clipping, omitted population,
hidden retry or shortened success. Independent zero-drift constant-sigma drag is
K0=(8/3)*sqrt(2/pi)*mu*sigma*sqrt(S); conserved energy/momentum also determine
eventual common temperature, without imposing it. Stated durations may precede
thermal equilibration, especially with physical electron/argon mass separation.

The outside-table probability is recorded throughout; it does not bound missing
force or heat, since the unsupplied cross section is unknown. The source's
stationary-target energy/angular convention and Maxwellian closure limits remain.
Outputs contain numerical evidence and scientific PNG/PDF, not experimental
validation or enabled moving-plasma neutral physics. Existing verification
infrastructure is reused; no unit suite is added or run.
