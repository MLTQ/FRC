# Directional accuracy, parent terms and strict arithmetic recovery

This continues the full experimental FRC sandbox goal. It remains a filled,
periodic, fully ionized verification fixture with physical electron/ion masses;
it does not select gas, fill pressure, coils or a formation experiment. Read
DIRECTIONAL_TERM_DESIGN.md for equations and contracts. The independent fixed
radial16,N128 reference and all older source/array archives remain unchanged.

## Separate spatial controls

The original balanced solver fully completed and replayed the axial256 and
fluid64 signed pairs. Physical relative volume L2 errors at quadrature20:

| Fluid radial / axial | Field radial / axial | vr rate | Br rate | J rate | Btheta rate |
|---|---|---:|---:|---:|---:|
| 32 / 128 | 128 / 128 | .108836499 | .051937972 | .087568776 | .383069049 |
| 32 / 256 | 128 / 256 | .108828829 | .002319415 | .084734261 | .259509336 |
| 64 / 128 | 128 / 128 | .112012914 | .051938079 | .084677158 | .380292248 |
| 32 / 128 | 256 / 128, split inverse | .030908386 | .051624975 | .033750429 | .279381110 |

The fluid64 pair gives vr .112012914, Br .051938079 and J .084677158: fluid-only
refinement does not resolve the field error. Axial refinement greatly improves
Br, but not the radial-momentum or current rate. Every signed input, coefficient,
physical fingerprint and non-runtime metric in these studies replays bitwise.

The original F32/128,B256/128 solve rejected the minus state. Its final first
attempt had operator2.604465e-10 and dual Faraday2.827549e-10, versus the original
1e-10 bound; physical Faraday was8.708597e-11. All five outer attempts failed.
The failed original case and its full bitwise replay remain failed records.
Exact input SHA256:
`77e43521a93bb9e51035655a1b3b7d1cdca26ab9c56fb6274c92d56603b32330`.

## Independently differentiated original terms

The new independent reference separates original ideal, Hall, pressure and
finite-current-inertia electric fields and their physical cylindrical curls,
plus pressure and Lorentz acceleration. Off-grid finite differences check each
curl, while summed fields check the independent mode equations. The production
capture evaluates the actual archived signed states, and its original fluid
terms reconstruct the complete saved rate bitwise. Original full Ohm remains
within the unchanged bound; no full rate or state is altered by the audit.

At F32/128,B128/128,epsilon1e-4, pressure and Lorentz radial acceleration norms
are421.787 and422.469 times the small net radial acceleration. Their own-term
relative errors are.000100923 and.000260343. Errors divided by the net rate
are.0425680 and.109987; complete vr error is.1088365. The outermost fluid element
carries92.50% of Lorentz squared error and89.11% of full radial-momentum squared
error. The explicit wall-reaction error is only.00675418 of net rate, so the
wall reaction alone does not explain the large error.

Toroidal induction illustrates stronger cancellation:

| Term | Analytic term norm / net rate | Term error / net rate | Own-term relative error |
|---|---:|---:|---:|
| Ideal | 144.227875 | .180962140 | .001254696 |
| Hall | 24324.080221 | 5.802628611 | .000238555 |
| Full pressure | 24408.196721 | 4.013111000 | .000164417 |
| Current inertia | .614170147 | 4.698210027 | 7.649688039 |
| Accepted net | 1 | .383069049 | .383069049 |

Inertia uses the numerical accepted J rate and is response feedback, not an
independent driving error. Error cross integrals are saved: Hall/reduced-pressure and
Hall/inertia errors correlate approximately-.548 and-.742. Summing individual
error magnitudes would be misleading. Only28.82% of full Btheta squared error
lies in the outer fluid element. The compatible-pressure-potential volume curl
noise is.000749963 of the net Btheta rate. Strong volume curl and explicit
wall-projected curl are separately recorded; volume curls omit face distributions.

## Driving-only control on the failed case

A separate capture evaluates original fluid and electric forcing at both valid
signed inputs, without calling induction or manufacturing a missing inertia
rate. Its shared baseline coefficients match the earlier full capture bitwise.
Every record explicitly retains the original full-solve status and accepts no
complete plasma rate. Quadrature20 errors normalized by the full reference rate:

| Fluid / field radial cells, axial128 | vr momentum | Hall toroidal driving | Pressure toroidal driving | Combined toroidal driving |
|---|---:|---:|---:|---:|
| 32 / 128 | .108836499 | 5.802628611 | 4.013110295 | 4.919187695 |
| 32 / 256, original solve failed | .030908386 | 2.232811195 | 3.090009085 | 2.078152578 |
| 64 / 128 | .112012914 | 5.802628611 | 3.666143267 | 4.683651591 |

The combined driving target excludes inertia. It is not an achieved magnetic
rate error. These measurements distinguish fluid and field contributions even
when the old full solve cannot be accepted. Both complete captures and physical
audits replay bitwise. Order16/20 physical volume term norms differ at most
7.65e-11 relatively in these six audited cases; all absolute scales, axis/wall
extrema, raw per-cell integrals and error cross terms remain in the records.

A further frozen reduction holds the physical wall band fixed at
31/32<=r/R<=1 across all meshes. Its fraction of full radial-momentum squared
error changes from89.11% at F32/B128 to2.263% at F32/B256, but increases to96.60%
at F64/B128. The corresponding Lorentz fractions are92.50%,3.954% and96.60%.
Thus the wall-dominated error is strongly reduced by field refinement. The
remaining F32/B256 pressure/Lorentz radial errors are.0425680/.0237568 of the
net rate; interior cancellation and axial accuracy still need work. These are
fractions of global squared error, not relative norms with a shrinking local
denominator. Fixed-band integrals and full error correlations replay bitwise.
Fluid-only doubling does improve the pressure radial error from.0425680 to
.00135002, a factor31.53, but leaves Lorentz error.112012 and therefore does not
improve the full momentum rate. This is evidence for separate fluid and field
requirements, not evidence that fluid refinement is ineffective in every term.

## Arithmetic diagnosis and recovered full input

Observation-only instrumentation reproduces the original first-attempt
rejection exactly. The proposed HT multiplication-association explanation was
disproved: its predicted gap is about2.7e-13 versus the observed5.22e-10.

Independent matrix-action reconstruction instead finds mass-inverse
distributivity error9.50e-16 relative to E, propagated to5.20351e-10 in the dual
residual, versus5.21883e-10 observed. The full predicted vector leaves7.06%
unexplained; no exact closure is claimed. All proposed-coordinate vectors and
every action are saved and replayed bitwise. No solve is accepted by this audit.

The separate opt-in SplitInverseWallCore evaluates
E=M^-1 f + M^-1(I Jt), sharing the first inverse with the unchanged original
RHS, rather than M^-1(f+I Jt). The physical linear equation, actual density,
masses, spaces, original A, all acceptance norms and iteration budgets are
unchanged. The complete saved failed input passes at the first iterate:
operator4.55059e-12, dual Faraday5.27777e-11, physical Faraday1.63177e-11,
Ohm5.23195e-16 and Ampere2.72639e-16. Restored original full-species Ohm is
2.53165e-16. The actual inputs remain bitwise unchanged.

The independently contracted fluid/current/magnetic powers are
[-3.12802502821e-11,1.16075305680e-11,1.96727194645e-11]W. Their residual is
about-2.50e-19W, or3.99063e-9 of the exchanged absolute powers. This measured
power defect is separate from strict solver residuals. The restored full
Faraday rate-relative diagnostic is2.72200e-7, while the original full backward
check passes at3.49241e-17. The old unprojected-curl defect is.00162764; it is
not silently reinterpreted as the declared wall Galerkin equation.

Analytic cylindrical polynomial full B/J/E and original work controls pass at
nr6/12/24; failed iteration exhaustion is atomic. **470 full regression tests
pass in233.11s.** Frozen regression source closure:
`3e6a80b63ec95a5fb43db59db78d8fac83f69b4959a7dc09606fdd39c254fc84`.

## Complete signed arithmetic controls and checkpoint

The new same-grid baseline and radial256 pairs at epsilon1e-4 and1e-5 all
complete. Every saved array and non-runtime metric of all three cases replays
bitwise. Maximum unchanged solver residuals are1.44528e-11,5.30071e-11 and
5.49747e-11. The old/new baseline inputs and full fluid rates match bitwise;
centered B/J/E coefficient changes are4.39e-11/5.86e-10/2.18e-16 relatively.

At radial256, reducing amplitude tenfold gives vr .030908444, Br .051624975,
J .033750497 and Btheta .279381837. Centered coefficient differences are
2.81507e-5 fluid,1.31783e-4 magnetic,1.08516e-3 current and4.56187e-8 electric.
These Euclidean coefficient differences are not physical reconstruction norms
and are not discarded because the main physical error norms nearly agree.
Order16/20 rate/initial norm differences are at most1.73e-15 absolutely. The
largest relative quadrature change,1.43e-5, belongs to an initial axis-density
error only7.84e-11; its absolute change is1.12e-15. Full scales are retained.

All continuation workers are terminal and collected. Thirteen full bitwise
verification groups and42 new source archives pass;23 preceding radial16,
59 earlier physical-interface,61 matched and32 wall archives remain byte
unchanged. Twenty-two earlier calculation source closures still match.
No previously tracked Python source changed. The user's modified poster is
preserved. All research and issue state remain local; no remote is configured.

The visually inspected six-panel report is
`outputs/radial16_directional_report_v2/verification.png` and `.pdf`.
The first report and its frozen source are preserved; v2 moves the legend
clear of the bars and explicitly identifies the baseline in panels B-D.
The complete source/replay/test/worker inventory is
`outputs/radial16_directional_diagnostics/checkpoint.json`.

## Remaining acceptance

Next compare F32/256,B256/256 with F64/256,B256/256, holding the same radial16,
N128 independent reference fixed. Retain amplitude, quadrature and original
strict-residual controls. Local pressure,
Lorentz and current-inertia accuracy remain separate from arithmetic recovery.
Then repeat complete signed trajectories with genuine mesh/time/residual
convergence and local charge, enthalpy and Poynting face budgets.

The frozen guarded trajectory source has not adopted the new arithmetic.
Heated moving exact vacuum, transport/ionization and validity, physical coil/end
conditions, converged formation and experimental controls remain unfinished.
No numerical atmosphere, diffusion, filter, artificial mass or energy repair
was introduced. Closed flux and axisymmetric response do not establish 3-D
stability or particle confinement. The full goal and acceptance beads stay open.
