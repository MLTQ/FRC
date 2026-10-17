# Finite-current collisionless slab: completed first physical evolutions

The new planar kinetic component completes the original20ns Ar+/electron slab
with physical masses, initial n=1e14m^-3,Ti/Te=.2/2eV, zero drift and field,
occupied±10mm and remote ends±.2m. Initially the exterior is exactly empty.
It evolves current and parallel heat directly, without the existing isotropic
fluid pressure assumption or the uniform kinetic-heat fixture's zero-current
constraint. These are component experiments, not the provisional pressure-range
lab discharge or an FRC. No molecular reaction or magnetic geometry is enabled.

## Equations and numerical assumptions

Solve dt f_s + v dx f_s + (q_s E/m_s) dv f_s=0 and
dx E=sum(q_s integral f_s dv)/eps0, with E(left)=0. Initially factorized
Maxwellians and a straight Ex leave transverse velocities unchanged: their
per-particle energy remains T_s0*eV. Parallel pressure and central third-moment
heat flux evolve from the actual distribution. No collisions, physical
viscosity, isotropization, wall emission, sheath model or external circuit.

Uniform finite-volume phase cells use conservative minmod linear reconstruction
and Strang stream/accelerate/stream. This limiter introduces numerical diffusion,
especially at steep fronts. Its measured energy drift is numerical heating,
not physical collisional heat. Initial Gaussian velocity cell averages are
integrated without renormalization. Initial P0 velocity representation itself
raises total energy by.347222%,.0868056%,.0217014% for64/128/256 velocity cells;
these offsets are separately reported from evolution drift.

Velocity support is explicitly finite at±8 species thermal speeds, omitting
Gaussian probability1.24419e-15 and second moment8.20805e-14. Actual x/v outflow
is measured, never clipped. This is not a cutoff-free kinetic distribution.
Spatial/velocity cell resolution and actual dt change together in the main
family. Every actual stage checks its CFL; no retries or thermal repair.

## Independent free transport

Three full20ns zero-force fixtures are compared with independently integrated
Gaussian characteristics f=phi(s)*1(|x-vT*s*t|<a), including exact cell-average
overlap cuts. This mathematical limit has diagnostic charge/E but no electric
force; it is not a self-consistent charged plasma experiment.

| x×v cells | Electron density relative L2 | Electron first moment | Second moment | Third moment |
| --- | ---: | ---: | ---: | ---: |
|1280×64|.0247691|.0159372|.0156737|.0160947|
|2560×128|.00953627|.00448948|.00495849|.00423838|
|5120×256|.00332967|.00165704|.00127642|.00115022|

Free energy/number changes remain below3.27e-14/6.43e-14. The heavy front travels
only13.9um thermally, so78–312um cells do not resolve it: ion first-moment error
increases2.17%→4.86%→9.72% over that family. Do not hide this with total density.
Separate full20ns ion evolutions at15.625/7.8125/3.90625um show first-moment
errors16.75%→9.03%→4.44%, and third-moment errors13.89%→8.12%→4.03%.
They use an explicitly embedded±12mm free-transport subdomain; actual boundary
distributions and all outflows stay exactly zero, so padding to±.2m is exact.
Velocity128→256 at3.90625um changes first-moment error4.44%→4.47%; the sharp
front's spatial diffusion remains limiting. This does not refine the coupled run.

## Self-consistent full trajectories and sensitivity

| x×v cells | Steps | Maximum relative energy drift | Maximum relative number balance error |
| --- | ---: | ---: | ---: |
|1280×64|480|.0282663|1.03e-14|
|2560×128|940|.0108157|1.36e-14|
|5120×256|1880|.00405479|1.69e-14|

Total momentum-balance errors are at most6.50e-27kg/(m s) for the main family,
including electric boundary stress. Exact reflection symmetry makes this a
weak check of local force accuracy. It is not a substitute for field/current
convergence. Particle energy lost at x/v boundaries is included in the energy
ledger; charge removal can also change field potential energy. Measured sigma8
ion velocity loss is below.000368 particles/m² versus2e12 initially, so that
effect is far below the observed numerical drift in these particular runs.

Successive full-support electron current L2 differences are20.59% then15.84%;
heat differences11.08% then4.13%; E differences3.70% then1.55%. Current is not
converged. Ion density differences also grow. At2560×128, halving CFL changes
current.364%,heat.126%,E.0411%, much smaller than main phase-grid sensitivity.
Extending velocity support to10sigma at fixed ds=.125 changes current.146%,
heat.0508%,E.0165%. That wider domain also tightens the spatial timestep; this
is not a pure velocity-domain comparison. Its omitted probability is1.524e-23.

The completed2048-label isotropic-fluid endpoint differs from finest kinetic
current/E by102.21%/27.50% in kinetic-normalized L2. Both remain numerically
unconverged and solve different closures. These numbers show material closure
sensitivity, not which solver is correct or an error against an exact solution.
The kinetic distribution develops nonzero parallel heat and unequal parallel/
perpendicular pressure. The inspected scientific plot keeps these distinctions.

## Evidence, references and next action

Main source SHA6721a90428179197fc54e282643c820bb9f1d55b6a6157e966b00a17593d34d5.
All eight full two-species runs and four full ion-front runs are terminal.
Directories `outputs/vlasov_slab_*` retain complete initial/final phase states,
every accepted timestep/budget,21 profile times, frozen sources and exact
characteristic comparisons. `vlasov_slab_comparison` holds comparison source,
JSON and inspected PNG/PDF. No unit tests or existing unit suite were run.

The conservation-form Vlasov equations and role of particle-flux reconstruction
are discussed by [Filbet, Sonnendrücker and Bertrand (2001)](https://www.math.univ-toulouse.fr/~ffilbet/Papers/paper03.pdf).
This implementation uses the explicitly derived CFL-limited minmod flux in its
companion, not that paper's large-timestep PFC reconstruction.

Next improve phase-space resolution or transport where current peaks and the
heavy free front remain unresolved, measuring actual complete trajectories.
Keep the reaction-source failure diagnosis as the primary coupled-run next step.
Do not install this independent kinetic model as a validated FRC heat closure.
Magnetized nonuniform collisions/heat, gas-specific chemistry, walls and
converged formation are still open under the unchanged overarching goal.
