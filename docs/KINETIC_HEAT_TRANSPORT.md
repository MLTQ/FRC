# Kinetic heat transport: completed uniform-background evolutions

The paired collision pulse still exposes bulk nonlocal transport scales.
Rather than install a local conductivity or an empirical flux limiter in that
pulse, an independent kinetic component now evolves heat in a controlled
periodic problem. This advances frc-transport.4, which remains open for the
actual nonuniform moving plasma, finite drift and physical edge/end treatment.

## Equations and scope

`plasma_bench.kinetic_heat` expands a small electron perturbation about a uniform
Maxwellian in orthonormal Hermite/Laguerre velocity polynomials. It includes the
full linearized electron-electron Landau operator, including field-particle
terms, plus velocity-dependent electron-ion Lorentz scattering. Ions are
stationary and infinitely massive in this explicit limiting problem. The
quasineutral electric field enforces zero perturbed density/current; it is not
a prescription for removing current dynamics from the formation simulation.

The collision matrix integrates the Landau weak form directly as weighted
transverse-gradient products. Gaussian centre/relative coordinates remove the
relative-speed singularity analytically. There is no velocity cutoff, collision
floor, fitted relaxation, state projection, artificial viscosity or flux limiter.
Exact density, energy and self-momentum invariants vanish in the analytic
integrand. Multiplication by parallel velocity uses the Hermite recurrence.
Matrix exponentials evolve every retained moment through each complete interval;
velocity resolution and eigenvector conditioning are recorded.

Read plasma_bench/kinetic_heat.md for equations and normalization. The parent
operator is discussed by [Pfefferle et al.](https://arxiv.org/abs/1701.08037).
The zero-current thermal-wave setting and limitations of reduced nonlocal
closures are discussed by [Brodrick et al.](https://arxiv.org/abs/1704.08963).
Neither reference makes this linear periodic calculation an edge or formation
model. It contains no neutral physics, ion heat flow, finite-amplitude front,
magnetic-field curvature, plasma escape or sheath boundary.

## Completed results, 2026-10-16

Initial calculation89659 exited0; source
4af72f04aab6d4f08b560f34c46167d81ebd17231934cebe14182d8c25bcb2c3,
outputs/kinetic_heat_relaxation. It evolved degrees8/12/16 at four collisionalities
and collisionless degrees12/20/28. Stronger nonlocal cases retained visible
velocity-truncation effects, so a complete refinement92079 was run, exiting0
in122.283s: outputs/kinetic_heat_relaxation_refined, source
acb28fcfc9834a0be3d29ca5e96eed07747f6a1b255394ee12580a99275e16cc.
Both scientific PNGs were inspected and PDFs saved. No unit tests were run.

Refinement uses degrees20/24/28 (121/169/225 moments) and the same physical
intervals, with a more accurate independent collisionless reference. The
standalone code is local branch e275ffe57 and main commit469927bc8. Only these
four new files were brought into main; all live pulse source closures remain
unchanged and the paired moving-collision branch is still otherwise unmerged.

The independent local Lorentz solution is obtained by integrating Maxwellian
speed moments under zero current, giving kappa=64*sqrt(2/pi)/Z in units
n*vT²/nu0. For the full Z=1 operator, the commonly quoted Braginskii coefficient
3.16 is an approximate reference, not a fitted target; the normalization used
here is documented in [Kupfer, Harvey et al.](https://compxco.com/Kupfer_Harvey_etal_1996_PhysPlasmas_3_3644.pdf).

| Velocity degree | Full Z=1 kappa / (n vT² tau_e) | Relative difference from analytic Lorentz kappa |
|---|---:|---:|
|8|3.1639876|-3.3207e-3|
|12|3.1881263|-1.9516e-4|
|16|3.1966766|-3.8686e-5|
|20|3.1998347|-1.2325e-5|
|24|3.2012848|-5.0780e-6|
|28|3.2020337|-2.4604e-6|

The full conductivity is increasing with velocity resolution; its degree28
value is1.33% above the approximate3.16 reference. No coefficient was rescaled.
For collisionless evolution over0<=k*vT*t<=4, an independent Gaussian-
characteristic Volterra solution gives maximum temperature differences
8.2171e-2,1.0610e-5,3.4723e-8 for degrees12/20/28. Halving the reference
integration interval from8000 to16000 intervals changes that reference by
1.0410e-7. The finest kinetic comparison is near the reference integration
floor; degree12 clearly recurs before the end and is inadequate. Collisionless
quadratic free energy changes by at most2.89e-15 relative to initial.

Define K=k*vT/nu0, vT=sqrt(T/m), nu0=n*4*pi*(e²/(4*pi*eps0))²*lnLambda/(m²*vT³).
Temperature/heat amplitudes below are normalized by the initial small deltaT/T.
At K=.001 the run spans five approximate local decay times; the other runs
span k*vT*t=4. These normalized problems are not selected lab times or lengths.

| K | Degree28 final temperature amplitude | Local approximate final amplitude | Max temperature change, degree24 to28 |
|---|---:|---:|---:|
|.001|.00632258|.00673795|8.6157e-5|
|.03|.48224379|.38654251|1.2545e-4|
|.3|.15788875|7.4468e-5|6.4789e-4|
|1|.10255304|1.7374e-14|1.6960e-3|

The corresponding maximum heat-amplitude changes are2.7472e-6,7.5326e-5,
5.5594e-4,3.1062e-3. Nonlocal relaxation is materially slower than local
diffusion. The strongest case still has about1.65% sensitivity relative to its
remaining final temperature amplitude; do not call all kinetic regimes fully
converged. Maximum eigenvector condition is8377 and imaginary temperature
roundoff is below3.25e-12. All sampled collisional quadratic free-energy
increments are negative. Periodicity fixes the spatial mean energy, so this
does not demonstrate a physical boundary heat budget.

## Next physical work

Use these complete kinetic histories as references for spatial heat evolution,
not as a local coefficient fit. A pulse-compatible model must retain history,
nonuniform thermodynamic profiles and changing field direction, finite current
and coupled ion dynamics. Explicitly derive how heat and any escaping particles
cross the moving plasma edge and ends; a periodic zero-current constraint cannot
be silently inherited there. Extend actual transport only after its applicable
limits and those boundaries are stated. Keep the ongoing full-pulse time/space
comparisons; linear heat agreement does not resolve their current sensitivity.
