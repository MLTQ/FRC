# Matched continuum response: measured checkpoint

The full FRC experiment-sandbox goal remains active. This is bounded periodic,
fully ionized argon verification with physical masses, 0.24 m length, 0.035 m
radius, n0=1e20 m^-3, 5/1 eV species reference temperatures and 5 mT bias. It is
separate from the user's provisional one-metre tube and neutral-fill envelope.
Read MATCHED_MODE_SOLVER.md for equations, boundaries, proposal algorithms and
limits. Full regression: **446 tests passed in 234.86 seconds**.

## Separate field directions

Eight complete signed instantaneous checks replay all coefficient arrays,
independent physical fingerprints and non-runtime metrics bitwise, at quadrature
16/20. Source **7bb77514a800c5d6e6513c5f06701d55e51bc5081f2d869d746cb0d2a5f2f125**.
At fixed fluid 16 elements / 128 axial points:

| Field nr × nz | Radial B rate relative L2 | Current rate relative L2 |
|---|---:|---:|
| 32 × 32 | 0.547443527 | 0.480470928 |
| 32 × 64 | 0.0139749510 | 0.0109873934 |
| 32 × 128 | 0.000244258701 | 0.000875416785 |
| 64 × 128 | 0.000244245758 | 0.000219550989 |

Radial 16/32/64 at fixed nz64 isolates the axial error floor. Additional controls
use epsilon 1e-4 versus1e-5 and eight versus sixteen fluid elements at 128 axial
points. All original physical rates, initial representation, axis/wall/max and
strict induction/power diagnostics remain. These rates alone do not verify a
complete trajectory on the fine meshes.

## Complete physical response

Six complete 100-microsecond trajectories form three opposite-sign pairs on
fluid elements4 / points32, field nr8 / nz32: 30 and 60 time steps at epsilon 1e-4,
and 30 steps at epsilon 2e-4. Source
**f255a57614954fdb125fde91e15fb4bf8e35abd8157d7e51a7f94255f2cd5b02**.
Every actual conservative state and Gauss stage is saved and replays bitwise,
along with every non-runtime iteration metric across all six trajectories.
All three independent
physical audits and both time/amplitude comparisons replay their full metrics
and physical fingerprints bitwise. The checkpoint records all six complete trajectory replays separately.

Selected final centered-response errors at quadrature20:

| Physical quantity | 30 steps | 60 steps |
|---|---:|---:|
| Density | 0.000968304 | 0.001010704 |
| Radial velocity | 0.002153991 | 0.002039019 |
| Radial B | 0.001907072 | 0.001936514 |
| Current | 0.002742451 | 0.002747908 |
| Full electric field | 0.008483819 | 0.008678388 |
| Electron pressure | 0.002998330 | 0.003150097 |
| Ion pressure | 0.003020999 | 0.002994500 |

The independent expected growth factor is 2.402175837. At 60 steps, physical
component projections range from 2.401218278 to 2.402203516. These are amplitude
diagnostics, not a fitted production eigenmode or a laboratory growth prediction.
Initial radial-B/current/E representation errors are 0.001105287/0.000586181/
0.006485440 and are reported separately. Endpoint errors need not decrease
monotonically when this spatial approximation dominates.

Halving the time step changes radial B/current/full E by 0.000282029/
0.000284198/0.001361556 relative to the first physical response. Two time
resolutions do **not** establish an order. Doubling amplitude changes these by
1.018715e-6/2.069890e-6/1.136489e-6. The amplitude effect is much smaller than
current coarse errors; full spatial and further time/residual controls remain.

Whole-domain L2 norm values at quadrature 16/20 differ by at most 1.0091e-15 in
absolute dimensionless value. All six states have positive independently sampled
density and temperatures. Stable original continuous-energy changes divided by
initial *total background* energy are at most 3.82255e-17 in magnitude; particles
and both entropy changes are retained separately. This normalization is not a
certificate of tiny perturbation-energy accuracy or local interface balance.
No field, heat, energy, physical mass or accepted state is repaired.

## Solver controls and retained failures

Uniform modal blocks independently match full physical blocks at rest/motion,
including every Fourier mode and Nyquist. Actual uniform nonlinear centered
RHS differences converge quadratically to the separate modal parent operator.
Tests independently substitute actual stages/endpoints and check transactional
failure, including preserved continuation history. Uniform current-inertia
preconditioning passes original physical-operator substitution and changing-n
RHS equivalence. All original strict induction checks remain.

The fine performance fixture uses fluid elements8 / points64, field nr16 / nz64.
Each pilot's accepted saved stages/endpoints pass a fresh full-equation audit,
and those audit metrics/fingerprints replay bitwise. Whole pilot trajectories
are performance evidence, not independently replayed convergence experiments.
RHS counts for successive accepted steps:

| Proposal variant | Step 1 | Step 2 | Step 3 |
|---|---:|---:|---:|
| Accepted-iterate secants | 191 | 163 | — |
| Continued secants | 191 | 214 | 40 |
| Continued Newton products | 143 | 512 | 430 |
| Bounded approximate proposals | 95 | 110 | 114 |

The unguarded coarse positive trajectory finishes, but its negative partner
rejects step 2 at residual 1.05044814e-12 versus original limit 1.00030777e-12 after
60 iterations. Its terminal failure and last accepted 3.3333-microsecond state
are preserved. Guarded proposals complete both signs without weakening this
bound. Earlier unarchived exploratory extra-tight atol 1e-14 work also rejected;
no passing 1e-14 control or bitwise evidence is claimed for that prototype.

One postprocessing worker had imported the old five-return fixture unpacker
before its correction and failed after the successful sixty-step runs. Its
failure directory explicitly records that the on-disk archive did not identify
the stale loaded module. No audit result was produced there. A fresh corrected
process produced the complete audited/replayed results above. Original simulation
source and trajectories were unaffected. The first plot layout is retained;
the corrected scientific PNG/PDF has been visually inspected.

## Evidence and remaining gate

Raw rates: outputs/matched_mode_directional_rates and matched_mode_fine_controls.
Trajectories/audits: outputs/matched_mode_guarded_ne4_z32_*.
Solver pilots: outputs/matched_step_*_pilot and their independently replayed audits.
Report: outputs/matched_mode_report/convergence.png and convergence.pdf.
Source preservation, regression and final replay status: outputs/matched_mode_diagnostics.
All 22 earlier verified calculation closures and all 32 preceding-checkpoint
source archives remain unchanged. No remote or publication is configured.

**frc-7mm.7 stays open.** Next: complete response mesh refinement (nr16/nz64 then
nr32/nz128 with matching fluid axial points), further time/residual controls and
physical small-signal energy scaling. Resolve local current/enthalpy/Poynting
interface budgets before adopting heated moving exact vacuum, transport,
ionization, coil/end boundaries and converged FRC formation. Keep frc-goal open.
Closed flux, equilibrium or axisymmetric growth cannot establish 3-D stability
or particle confinement. The user's lab envelope remains an exploration input.
