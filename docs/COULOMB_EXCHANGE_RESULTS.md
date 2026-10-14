# Electron-ion energy exchange: complete plasma evolutions

The new opt-in CoulombExchangeWallPlasma adds conditional Maxwellian thermal
equilibration to the moving finite-mass solver. Equal and opposite heat sources
act on the actual electron/ion entropy variables. Physical masses, flow,
pressure/current coupling and magnetic evolution remain active. Neither total
energy nor temperatures are reset after a step. Both species require gamma=5/3;
the legacy contact fixture's gamma_i=1.4 is deliberately not reused.

The collision rate uses the SI Landau coefficient and the electron-ion
Coulomb-log branches documented in the
[NRL Plasma Formulary 2023, pp.33–34](https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf).
This closure assumes weak coupling, Maxwellian species, Z=1 and small relative
drift. Scalar resistivity and its heating split remain separately prescribed;
this is not a complete transport tensor, neutral or ionization model.

## Actual full-system runs

Uniform fully ionized argon: n=1e20 m^-3, Te=5 eV, Ti=1 eV, radius0.035 m,
periodic length0.24 m, zero flow/current/field and zero prescribed resistivity.
Both actual evolutions finish400 us. The independent reference separates the
scalar temperature ODE and solves its time integral along Te+Ti=6 eV; it
never calls the production collision rate or plasma RHS. Final reference
temperatures are3.000142079745 and2.999857920255 eV, approaching3 eV each.

| Measured result | 40 steps | 80 steps |
| --- | ---: | ---: |
| Final electron temperature error, eV | 2.05219e-6 | 1.30405e-7 |
| Final ion temperature error, eV | 2.03512e-6 | 1.29339e-7 |
| Maximum electron temperature error, eV | 3.61903e-6 | 2.28793e-7 |
| Maximum physical total-energy relative drift | 6.81217e-7 | 4.32892e-8 |
| Peak own RSS, bytes | 508444672 | 576602112 |

Halving dt reduces final errors by15.737/15.735 and energy drift by15.736,
consistent with the fourth-order Gauss method. These are two measured time
resolutions, not a spatial-convergence result. Independent physical sampling is
q20. Maximum spurious flow is1.07e-8 m/s, current1.68e-19 A/m² and field
2.69e-28 T across both runs. Entropy increases as expected. The saved particle
metric is a departure from the analytic uniform integral (at most2.64e-14),
including constant-field reconstruction error; it is not an isolated mass drift.

The initial SI energy-exchange frequency is6636.819411/s. The NRL rounded
cgs/eV prefactor1.8e-19 gives6809.497163/s, a2.536% difference. That coarse
coefficient comparison is separate from the high-accuracy scalar time integral;
it is not presented as a numerical temperature error or experimental validation.

Outputs are `outputs/coulomb_temperature_relaxation_t40` and `_t80`, including
every accepted state/stage, solver history and source archive. Actual sessions
20930 and96369 exited0 and were collected. The comparison and visually inspected
plot are in `outputs/coulomb_temperature_relaxation_comparison`. No unit tests,
new verification framework or automatic replays were added/run.

## Scope and next physical use

The pointwise collisional source conserves native thermal work and produces
nonnegative native physical entropy. Nonuniform reconstruction, strong drift,
heat conduction, physical friction/heating consistency, neutrals and ionization
remain to be evaluated. A next direct simulation is an advected two-temperature
contact with constant total pressure, compared with scalar relaxation along
its characteristics. Do not substitute these uniform results for moving-edge,
exact-vacuum, real coil/end or FRC-formation evidence.

The separate numerically damped contact trajectory has its own equations and
status in HYPERDIFFUSIVE_CONTACT_RESULTS.md. No new collision term was inserted
into that running worker or its frozen source closure.
