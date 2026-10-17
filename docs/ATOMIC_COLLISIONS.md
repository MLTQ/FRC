# Reacting charge populations with Coulomb current relaxation

## Moving planar coupling, 2026-10-17

`VectorizedCollisionalAtomicMoments` now combines the existing charged-pair
collision split with array-assembled reactions in the precise moving-support
solver. Updated populations, momenta and energies feed actual motion, pressure,
Gauss field and current. No duplicate prescribed resistance or heat is added.
The original homogeneous studies below retain their stated scope.

Run32422 is EXIT0 in `outputs/reactive_material_argon_collisional_heavy_10ns`:
full20ns,971steps (527/444),21349 final density cells,89.884s and Erel2.84193e-11.
It uses32 initial cells, variance return, heavy-only mechanical enrichment,
128-bit coordinate scaling,10ns source intervals and explicitly lnLambda10.
Source SHA `98867099cce773f723371386ea1095199ac4efeb6c667b31680b2a3bbeaf3e1d`.
Against the otherwise matched reaction-only scaled case, E/J relative differences
are1.56447e-6/8.96056e-4. Electron front shift8.56890e-10m, common-volume
velocity/T RMS0.71199m/s/0.000395534eV; these differences include moving sharp
features and source splitting, not a direct measurement of integrated drag heat.
Neutral T RMS9.039e-8eV and Ar1 T RMS8.752e-8eV are also recorded.

This is the low-density n=1e14m^-3 planar component with unresolved current and
thermal resolution. It does not establish collisional evolution of the chamber
fill or justify the chosen Coulomb logarithm. Electron/ground-neutral elastic
exchange is a separate explicit addition documented in
[ELASTIC_NEUTRALS.md](ELASTIC_NEUTRALS.md). Ion-neutral scattering, charge exchange,
excited-target scattering, screening and driven changing-species coupling remain.

## Earlier homogeneous source

The new local source evolves all retained charged-level and electron momenta and
energies with elastic Coulomb exchange and the supplied atomic reactions. The
equations, conservative pair construction and explicit limits are documented in
plasma_bench/atomic_collisions.md. It uses the same isotropic finite-drift Landau
moments already compared with independent velocity integrals; lnΛ=10 is a stated
input, not an inferred screening model. Neutral collisions, charge exchange,
spatial fluxes and electromagnetic feedback remain absent.

## Preserved first attempt

Run 11914 exited 1 in 1.969 s, source
73c68c8d711f4a8c7f5d7e807b03c763072c3f766de170b8213f6ae99db2bbc3.
All three collision-only histories completed and converged toward the archived
independent velocity-space reference. Both reacting 5/2.5 ns cases rejected their
first step. The actual rejected 2.5 ns state is preserved with the full source
in outputs/atomic_collisions. Charge 9..15 populations (3.0e-81 down to 3.3e-294)
had negative thermal energies. This is not a subnormal-only rounding problem:
the frozen collision energy representation was advanced across creation of new
populations. No tiny population was removed and no energy was overwritten.

The revised step is collision-half / reaction-full / collision-half. Populations
remain constant during each collision substep, and both moments and energy from
each source enter the other source. Every substep still validates actual thermal
positivity; the method is second order, not unconditionally positive. All physical
reaction/collision equations and parameters remain unchanged.

## Independent collision limit

The existing electron/argon velocity-space reference is reused at all archived
times. It integrates binary transfer in velocity space without the production
erf/incomplete-gamma formula. The matched problem lasts 1.088042847 µs, begins at
three combined thermal speeds of relative drift, and uses physical electrons
and the reference's 39.948 u ion mass. Atomic conversions are explicitly disabled
in this limit. Its high drift does not validate the separate small-drift reactions.

| Full outer steps | Maximum drift error / initial drift | Maximum electron-temperature error, eV |
|---:|---:|---:|
| 100 | 2.33708e-4 | 0.0103561 |
| 200 | 5.94782e-5 | 0.00284754 |
| 400 | 1.53257e-5 | 0.000724248 |

The decreases are approximately second order. Finest ion-temperature error is
1.48862e-8 eV. Maximum total energy-budget drift across these complete limits is
7.95808e-13 eV/nucleus and scaled total momentum drift is <=1.67833e-14. These
quantify the collision integration, not complete multispecies/field accuracy.

## Complete heated argon source

Run 19857 exited 0 in 311.275 s in outputs/atomic_collisions_split, source
acc4681d50093851137f8254b2f20528b69e76f6a4439a711d21db0137b3f782.
Both 700/1400-step heated argon histories complete. They use the previous physical
3.5 µs protocol: 10% initial
seed, Te=10 eV, heavy T=0.1 eV, N=8.961429541213399e19 m^-3 and 20 eV/nucleus
prescribed heating. The neutral mass is 39.948 u with charge-dependent electron
mass subtraction. Initial neutral/ion axial velocities are +1000/-1000 m/s;
electron azimuthal velocity is 10000 m/s. Fixed frame is the laboratory frame.

The finer trajectory ends at 0.683265616 electrons/nucleus, Te=2.502600597 eV,
number-weighted ion T=0.263221409 eV and neutral T=0.100004771 eV. Current
magnitude falls from 14429.4034 to 0.016713821 A/m²; the final vector is
[0,-1.14761e-6,0.016713821] A/m². It first reaches half and one tenth of the
initial magnitude at recorded times 70 and 200 ns, and one hundredth at 472.5 ns.
The reaction-only control retains roughly 14429 A/m² and ends at Te=2.536280756
eV. This is collision-driven local relaxation with no sustaining electromagnetic
drive; these are not decay predictions for the driven laboratory plasma.

Time halving changes population L1 by at most 1.72795e-5, electron T by
0.001049757 eV and the current vector by 0.616047 A/m² (4.27e-5 of initial
current). The small final residual current itself differs by about 3.74%; a small
error relative to initial current does not establish that tail accurately.
Maximum energy-budget drift is 4.50795e-12 / 2.54685e-12 eV/nucleus in the
coarse/fine histories; maximum nucleus/mass drift is <=2.724e-13, charge error
<=9.304e-14 and scaled momentum drift <=1.627e-13. The small relative-drift
reaction ordering reaches 0.0053584 here. These are homogeneous numerical
measurements, not experimental validation or whole-sandbox convergence.

One coarse accepted state at 410 ns has an Ar17+ population of 5e-324/nucleus
and represented zero thermal energy; it exercises the stated one-cold-species
limit without a floor. The finer run has no cold occupied charged species.
This finite-representation observation is recorded in physical_readout.json;
it is not evidence of physically cold material in the experiment. The full
scientific PNG was inspected and its PDF is preserved with both histories.

All accepted states and original failures are retained. Raw personal provider
tables and rate matrices are not bundled. All missing atomic channels, branching,
photon transport, distribution validity and wall/edge issues remain open. This
source is not yet attached to the fixed-carrier moving-coil state; that requires
actual species transport and changing current inertia, rather than resetting n,
temperature or current after an atomic update. The complete FRC goal stays open.
