# Complete weak entropy budget and fixed-current refinement — 2026-10-12

This continuation from 1146f0b made numerical progress. The full research goal
and all acceptance beads remain open. No physical solver, accepted state,
mass, tolerance or boundary condition was changed. No unit tests were run.
The inspected figure is [species_budget_report/verification.png](../outputs/species_budget_report/verification.png).

## Whole-vector interface identity

The new independent root-polynomial/Fourier tensor adjoints reconstruct every
coefficient of the original species K action. The decomposition retains volume
advection, half the physical current divergence, both one-sided face fluxes,
and the continuous-potential transfer residual. Axis/wall and periodic seams
are evaluated separately. It is an integration-by-parts identity, not a new
numerical interface flux. Both species and the saved-error/independent-mode
directions are checked at quadrature 12, 16 and 20, with original and independent
potential projections separately retained.

Maximum full-vector sum-versus-direct relative difference is 1.78627e-12.
The independent direct action agrees with the original production action
within 2.183e-12. Independent q16/20 currents, divergence and every one-sided
trace reproduce the preceding charge audit exactly. Independently integrated
potential maps agree with production maps within 2.649e-14 radial and 5.346e-15
axial. No quadrature-induced change of the operator is concealed.

The analytic control Jr=r*k*cos(kz), Jz=-2*sin(kz), divJ=0 and
s=r²/R²+0.2*cos(kz) verifies the full-vector identity within 4.775e-14 and its
analytic advective covector within 2.162e-13. This control permits nonzero normal
wall current; it is separate from the physical wall fixture.

For the saved-error electron entropy at q20, the physical directional
contributions (/s) are:

| Contribution | Projection |
|---|---:|
| Volume advection | +950243.809990 |
| Half volume divergence | -994306.301502 |
| Radial interior faces | +1143532.174065 |
| Potential transfer residual | -149018.773177 |
| Axial interior faces | -206.855592 |
| Periodic seam | -0.243756 |
| Complete sum | +950243.810027 |

Axis contribution is zero and the wall contribution is near floating roundoff.
These are signed projections, not eigenvalues, causal fractions or additive
error norms. Large terms cancel. The complete physical rate differs from the
volume-advection covector divided by the actual fluid point volume by only
2.03417e-8 relative for electrons and 1.74350e-8 for ions. In the smooth control,
those differences are 1.90652e-7 and 1.61521e-7.

The archived centered electron rates approach the exact stationary linear
expression as amplitude falls: saved-direction physical differences are
7.47736e-9, 7.47744e-11 and 7.49413e-13 at h=1e-4,1e-5,1e-6. Both signs and
both species are retained. This does not certify nonlinear evolution.

The base representations S(C/n) and (S C)/(S n) differ by only about 3.3e-14
relative; axial derivatives agree within 4.3e-14. Their nominally zero radial
derivatives give ill-conditioned relative comparisons of hundreds, with maximum
absolute differences only 4.147e-10 per meter. Those ratios are not evidence
of a large physical radial gradient. The representation distinction therefore
does not explain the observed 55% strong-source discrepancy in this base.

## Physical mass and source representation

The next audit independently integrates M=S^T W S and retains the original
diagonal volume V. It decomposes the actual physical discrepancy into:
production-minus-volume/diagonal rate; diagonal-minus-exact-mass rate; and
exact-mass projection-minus-strong-source. Their full fields and pairwise
signed products are retained; their separate norms must not be added.

| F8 saved-direction relative norm | Electron | Ion |
|---|---:|---:|
| Complete error against strong source | 0.550011544 | 0.549686500 |
| Interface/transfer remainder | 1.66486e-8 | 1.42470e-8 |
| Diagonal-versus-physical mass effect | 0.128754250 | 0.134411016 |
| Source representation remainder | 0.534728942 | 0.532999932 |

Replacing only the mass in this diagnostic barely reduces the full error;
most of the source is unresolved by the original fluid space. No mass change
is adopted. The independently integrated normalized mass matrices agree with
analytic Legendre/Fourier projectors within 4.183e-14 radial and 3.219e-15 axial.
The highest local degree-five radial mode has mass ratio 5/11; the real axial
Nyquist mode has ratio 1/2. Their derivation is recorded in
SPECIES_MASS_DIAGNOSTIC_DESIGN.md. The saved source has 15.83%/16.98% of its
point-volume-normalized coefficient energy in the highest radial subspace for
electrons/ions; this metric is distinct from physical entropy-rate energy.

Native skew/constant conservation products are near roundoff when normalized
by the corresponding absolute products. Physical entropy-density-rate integrals
divided by their absolute integrals stay below 9.109e-14 across the mass audit.
Unnormalized residuals can look large because these direction amplitudes are
arbitrary; they are not laboratory entropy loss rates. This is instantaneous
global accounting, not the unfinished local/time-integrated energy budget.

## Refine the fluid space while preserving the current

The B32/64 current coefficients and their physical reconstruction are unchanged
in every case. The original stationary base is bitwise constant in radius and
is exactly replicated onto each new fluid grid. All field cells still subdivide
fluid elements. The F8 original, diagonal-advection and exact-mass coefficients
reproduce the previous audit bitwise. No accepted state is transferred/evolved.

Errors against the same original strong entropy source at q20:

| Fluid radial elements | Electron, original rate | Ion, original rate | Electron, exact-mass diagnostic | Ion, exact-mass diagnostic |
|---:|---:|---:|---:|---:|
| 8 | 0.550011544 | 0.549686500 | 0.534728942 | 0.532999932 |
| 16 | 0.089749439 | 0.091143079 | 0.088836939 | 0.090154604 |
| 32 | 0.000292433 | 0.000294120 | 0.000011637 | 0.000009848 |

The smooth independent mode stays near 2.6434e-5 electron and 2.3402e-5 ion
under radial refinement; axial and field spaces remain fixed. Actual floating
reconstruction changes in the nominally unchanged base source are reported,
at most 7.554e-12 relative, rather than hidden as identical sampled data.
This is scalar-source spatial convergence, not full coupled mesh convergence,
an eigenvalue study or proof that the growing trajectory is cured.

## Reproduction and resource evidence

All actual worker exits were collected as zero. Full fresh-process replays
match every raw byte and all non-runtime metrics: 9 weak-budget files, 4 mass
files, 12 space-refinement files and the report's integral file. The report PNG
also replays exactly. Its 128 q16/20 comparisons use actual field differences
on a common q20 grid; maximum difference relative to the full production rate
is 1.93384e-12. The largest component-relative difference is 1.79770e-11 for
a small axial-face term, whose full-rate-relative difference is 1.74067e-14.

The six guarded calculations/replays stay below 696647680 bytes peak RSS and
above 5939200000 bytes free disk, including final observations. Report stage
observations are recorded separately and are not a continuous guard. Read
`outputs/species_interface_audit`, `species_mass_audit`, `species_space_audit`,
their replay/verifier folders, and the corresponding diagnostics controls.
All preceding frozen sources and the user's original poster are preserved.

The final aggregate record is outputs/species_budget_diagnostics/checkpoint_v2.json,
with 200 Python sources. Its v1 record is retained: that initial aggregate seed
list omitted four dynamically dispatched scientific modules, although every
individual calculation archive already included them and was verified. The
v2 union is built from those complete individual closures; no calculation or
result needed rerunning for this aggregate-record correction.

Next: COUPLED_FLUID_REFINEMENT_DESIGN.md. Evaluate the complete unchanged
coupled RHS on exactly transferred directions before choosing a new evolution
or a consistent operator change. Full signed/time/arithmetic trajectories,
local energy, heated exact vacuum, transport/ionization, physical coils/ends,
converged formation and the experiment UI remain required.
