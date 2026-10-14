# Moving contact with numerical short-wavelength damping

The opt-in sixth-order damping candidate completed100 us, but does not reproduce
the independent continuum perturbation accurately enough for adoption. It delays
the earlier rapid growth; the error returns late in the run. The unresolved
moving pressure/current problem remains tracked in frc-7mm.7. No formation or
lab-design claim follows from solver convergence.

## Equations and independent target

HyperdiffusiveWallPlasma adds numerical magnetic/current damping and paired
numerical electron heat. Physical resistivity, ideal/Hall/inertial/pressure
terms and original solver acceptance bounds remain separate and unchanged.
Read the model, modal-solver and contact-example companions for equations,
variational higher-order wall conditions and the numerical heat assumption.

The prospective mesh-based strength0.25 gives eta6=1.1876780167018035e-19
 ohm meter^5 at F8/64,B32/64. Electron numerical heat fraction is1. This is the
original filled periodic 70mm diameter,240mm long argon contact,5mT bias and
negative1e-4 initial perturbation, not the proposed laboratory chamber.

The independent Fourier/Bessel matrix exponential with the analytic added
eta6 curl^4 J term predicts changes at100us of1.81800e-7 for density,
6.61795e-7 for Btheta and4.81299e-6 for Bz, relative to each continuum signal.
The unchanged exponential agrees with the original growing mode within3.91e-11
in every component. Raw operators/vectors are in
outputs/hyperdiffusive_contact_reference_s025. This predicts resolved-mode
continuum distortion, not the modified discrete solver's accuracy.

## Complete nonlinear trajectory

outputs/hyperdiffusive_contact_f8_b32_s025 completed all60 steps through100us.
Actual session40062/PID41757 exited0 and was collected. The adjacent _control
record is terminal: elapsed2814.49s, peak ownRSS4800626688 bytes, minimum free
disk10280960000 bytes; no resource guard triggered. Every accepted state/stage
and source archive is retained. No scientific worker remains live.

The existing independent density observer covers every endpoint. At18.333us,
error was3.30690e-5 of the predicted perturbation versus0.0696199 in the otherwise
identical original undamped run. At31.667us it was6.42403e-5. It then grows to
7.54987e-4 at66.667us,0.157514 at91.667us and0.954513 at100us. The original
undamped run ended after28.333us with error790.258. These comparisons show a
substantial delay, not elimination of the accuracy failure.

Independent final physical reconstruction gives:

| Quantity | Error / corresponding continuum perturbation norm |
| --- | ---: |
| Density | 0.954513 |
| Radial velocity | 0.182582 |
| Electron entropy | 1.49660 |
| Btheta | 1.41193 |
| Bz | 0.541441 |
| Electron pressure | 2.00007 |
| Current vector | 25.9242 |

The reference is linear, not an exact nonlinear solution. The finite amplitude
and single mesh/time step remain limitations; these large discrepancies do not
justify adopting the candidate. All reported field ratios agree closely at
q16/20. The full density observer and endpoint data are in
outputs/hyperdiffusive_contact_density60 and hyperdiffusive_contact_endpoint60.
Observer81511 exited0; peak ownRSS1533788160 bytes. Earlier observers40/55
(sessions94697/71800) also exited0. No unit tests or automatic replays ran.

## Conservation and interpretation

At18.333us the independent physical energy drift was-1.23013e-7 of the initial
perturbation energy. At100us it is-7.87834e-13J, or-0.0476713 of that initial
1.65264e-11J scale (-5.01265e-12 relative to the full background energy).
Particle-relative drift is-1.24425e-17. Density and temperatures remain positive.
The quadratic perturbation energy is a normalization, not a conserved invariant.
Instantaneous native work pairing does not establish continuous physical-energy
accuracy or finite-step conservation. No stored state or energy was repaired.

The inspected plot in outputs/hyperdiffusive_contact_comparison displays the
whole density error history. It is a numerical comparison, not an experiment.

Metadata erratum: the original source/output scope calls radial cells nonuniform
in r. The actual field grid is uniform in r and nonuniform in x=r^2/2;
h=R/nr is the correct radial spacing used in the coefficient. Frozen source and
raw metadata remain preserved with this explicit wording correction.

## Next physical work

Do not repeat this same candidate or expand replay/reporting infrastructure.
Retain the unresolved pressure/current accuracy issue. The separately completed
Coulomb thermal-exchange component (COULOMB_EXCHANGE_RESULTS.md) can next be
used in a genuinely advected, nonuniform monatomic two-temperature contact,
compared with scalar relaxation along exact characteristics. That advances the
physical model without claiming that thermal exchange cures the magnetic/current
problem. Exact vacuum, physical coil/end coupling, consistent collision/friction,
conduction/ionization, converged formation and the full original goal remain open.
