# Moving resistive plasma — 14 October 2026

The complete moving finite-mass wall model now has an opt-in scalar resistivity
and a paired Joule-heating source. The old collisionless model is unchanged.
See `plasma_bench/resistive_wall_plasma.md` for the equations and limitations,
and `examples/resistive_bessel_decay.md` for the independent decay problem.
No unit tests or new replay/reporting framework were added.

## Completed end-to-end result

Two actual nonlinear runs evolve a small axial Bessel field for five analytic
decay times, 5.260160724 us, with all motion, species pressure, current inertia,
Hall and drift-stress terms present. Both completed all40 steps at the original
Gauss1e-12 absolute/1e-10 relative checks and original induction bounds.
Actual sessions44111 and18880 returned exit0.

The independent linear decay rate is950541.297537 /s, including physical
electron/ion reduced-mass inertia. The field has decayed to exp(-5) of its
original amplitude at the comparison endpoint. Finite-amplitude nonlinear
motion is present, so this linear target is an approximation, not an exact
nonlinear reference.

| Physical comparison | Field nr8, fluid4 | Field nr16, fluid8 |
|---|---:|---:|
| Final Bz L2 error / remaining analytic signal | 3.22082e-5 | 2.60917e-6 |
| Final Jtheta L2 error / remaining analytic signal | 5.01101e-6 | 4.62034e-6 |
| Maximum native total-energy drift / initial field+current energy | 4.64017e-11 | 9.28034e-11 |
| Maximum independent physical total-energy drift, same denominator | 4.64017e-11 | 9.28034e-11 |
| Peak physical speed, m/s | 0.0191226 | 0.0191228 |
| Peak own RSS, bytes | 434290688 | 847380480 |

Independent root-polynomial/Fourier quadrature20 reconstructs the actual
fields and EOS. The initial analytic magnetic-plus-current energy is
5.981582722e-7 J. At nr16 the thermal-energy gain is0.999954334 of that energy,
bulk kinetic energy is2.21863636e-7 of it, and remaining magnetic/current energy
is4.54001136e-5 of it. Initial discretization error is retained separately.
The tiny total-energy differences approach floating subtraction resolution
against the large thermal background; a printed zero is not exact conservation.
Mesh refinement improves B strongly, while current error changes much less.
Time, finite-amplitude and initial representation effects are not separated
by these two runs; no formal convergence order is asserted.

The model uses explicitly prescribed eta=1e-4 ohm meter and electron-only
deposition for this experiment. It does not yet calculate a temperature- or
gas-dependent resistivity, temperature exchange, ionization, radiation or heat
conduction. High-order projected heat is not guaranteed pointwise positive
for unresolved profiles. The experiment's fixed cylinder and gamma_i=1.4 are
numerical fixture choices, not laboratory parameter recommendations.

Evidence: `outputs/resistive_bessel_decay_nr8` and `_nr16`, each containing
the complete initial/accepted states, stages, source archive and measured
history. Physical positivity passed every stage and endpoint; no floors,
extra numerical diffusion, changed masses or energy corrections were used.

## Collisionless limitation remains

The F32/64,B32/64 quarter-step contact run reached accepted states but failed
the intended physical comparison: at13.333us, density error was0.0482122 of
the predicted perturbation versus0.000659527 for the earlier F8/64 trajectory.
Early five-part-per-million agreement did not persist. Fluid resolution and
time step both differ, so this is not isolated spatial/time convergence.

That run was deliberately interrupted after observing the growing discrepancy;
all accepted states remain. Read `scientific_stop.json` in
`outputs/trajectory_fluid32_t240_minus` and COUPLED_FLUID_REFINEMENT_RESULTS.md.
No100us pass, stage-rejection result or collisionless cure is claimed.
The successful resistive decay problem does not excuse this unresolved mode.

Next simulation work should investigate the unresolved pressure/current coupling
and extend the new resistive dynamics to
nonuniform moving cases, using existing independent continuum comparisons.
Exact vacuum, physical coils/ends, transport/ionization and formation remain
open. Do not make replaying these completed runs the default next milestone.
