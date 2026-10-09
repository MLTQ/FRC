# Fixed radial16 reference, spatial response and strict solve recovery

Measured 2026-10-09. Design, equations and scope are in
[RADIAL_MODE_DESIGN.md](RADIAL_MODE_DESIGN.md). The full laboratory sandbox goal
remains open. These are filled-periodic instantaneous verification fixtures.

## Independent reference refinement

At one fixed radial wavelength, the independent growth rates are
12307.156241463335/s (N64), 12307.15614451037/s (N128), and
12307.156144548411/s (N256). Recomputed N64 raw arrays exactly recover the old
`stationary_contact_audit/continuum_07.npz` SHA256
`f9c92d8d94be18b10955fdde1f20868e2d3f0fafcc58a7157eca3e318258f6c6`.

Growth convergence hid a physical-reference limitation: the original even-grid
Fourier convention retains a real Nyquist cosine in magnetic values but sets
its collocation derivative to zero. Off-grid fourth-order finite differences
find Jr curl mismatch .00325101825 at N64, unchanged over three difference
steps. N128/N256 reduce it to about1.3e-10/8.1e-11. After diagnostic scalar
alignment, N64-versus-N256 electric-theta shape error is .00588697035;
N128-versus-N256 physical/J/E component errors are all below6.9e-9. Independent
quadrature12/16 and full array/metric replays are retained. No coefficient was
filtered and no old source changed.

All following controls use the same N128 reference, raw SHA256
`3b142fe986aca23d340b1739c17ec83b5b94e04574cff2a661fe1b0a8527ef08`.
Outputs: `radial16_reference*`, `radial16_reference_audit*`.

## Separate fluid and field resolution

F is fluid radial elements / axial points. B is field radial cells / axial
cells. Entries are relative physical L2 errors of the full centered signed RHS,
at quadrature20 and amplitude1e-4. They are fractions, not percentages.

| F | B | vr rate | Br rate | J rate |
|---|---|---:|---:|---:|
|4/64|16/64|344.215|151.221|476.613|
|4/64|32/64|189.911|119.738|792.385|
|8/64|32/64|20.6316|8.29375|71.1844|
|8/64|64/64|15.7960|2.80474|65.6822|
|16/64|64/64|3.86960|1.52766|5.20794|
|16/128|64/128|3.86912|.129238|5.10160|
|16/128|128/128|.956840|.0547405|1.01045|
|32/128|128/128, balanced solve|.1088365|.0519380|.0875688|

The finest complete case also has density .00243961, vtheta .000137404,
vz .0000428128, electron entropy .000398270, ion entropy .000103647,
Btheta **.383069**, Bz .00549650 and full E .000140880. In particular, a
38.3% toroidal-field-rate error remains. Refinement has not closed the accuracy
gate. No convergence order is inferred from this mixed refinement table.
The finest wall-region vr/Btheta/Jt errors are .576775/1.154478/.438075;
axis-region values are .001189/.274473/.032686. Wall-relative norms use the
reference norm in that region; they are not fractions of the global error.

Initial representation converges much faster: final density/vr errors are
1.05346e-5/1.04147e-5, Br/Bz 7.11922e-5/.000116518, weak J .000266164.
The strong-volume curl(B)/mu0 versus weak-current discrepancy is .000982571,
normalized by the **strong** current norm; it excludes surface distributions.
This differs from the denominator of the previous extracted-mode .599823 metric.

Radial normal-current jumps / analytic face norm decrease
3.201116 -> .197996 -> .00919278 -> .000320212 on16/32/64/128 field cells.
Radial tangential-B jumps similarly reach .000323249. At128 axial cells, axial
normal-J/tangential-B jumps remain9.93422e-5/.000126339. Both physical sides,
per-face integrals, max/axis/wall errors, hashes and quadrature controls are saved.
Outputs: `radial16_matched_rates*`, `radial16_fine_rates*`,
`radial16_balanced_rates*`. All failed records remain explicit failures.

## Amplitude, quadrature and baseline controls

At the finest grid, epsilon1e-5 gives vr/Br/Jt rate errors
.1088361693/.05193797165/.08756906320, compared with
.1088364987/.05193797164/.08756877589 at1e-4. The errors are nearly unchanged,
but the centered coefficient vectors are not identical: relative Euclidean
fluid/B/Jt/E differences are2.81486e-5/.000115461/.000773103/2.13268e-8.
These coefficient norms are explicitly distinct from reconstructed physical
norms and are not a proof of fully resolved linear amplitude response.

The largest quadrature16/20 relative discrepancy in a volume L2 norm is
3.04670e-11. Across axis/wall norms the relative maximum is1.42756e-5 for an
initial axis density error of only7.84e-11; its absolute difference is1.12e-15.
The largest absolute norm difference is4.55e-12 in the very large coarse wall
Btheta-rate error. These numerical scales are recorded rather than concealed
by a single universal relative tolerance.

The already accepted F16/128,B64/128 baseline matches **every saved array
bitwise** between old and new induction variants. The amplitude/quadrature
comparison itself replays all arrays and metrics bitwise in
`radial16_comparison*`. The PNG/PDF `radial16_report/verification` has been
visually inspected and records its exact input-summary hashes.

## Preserved strict failure and opt-in recovery

The original F32/128,B128/128 minus solve exhausts every outer attempt; its
operator residual remains1.39069e-10 under the original1e-10 limit, while
physical Faraday is2.62386e-11. The exact failed input archive SHA256 is
`e7a3cf637ceff0f3201a88f12ddd84b593f28d23f898406915c8d040a3eac755`.
Its complete error and bitwise failed replay remain in `radial16_fine_rates*`.
No rejected rate was used in a measured passing case.

New opt-in balanced refinement preserves all physical equations, masses,
preconditioning, three inner/five outer budgets and every acceptance bound.
On the exact captured input its original operator / dual-Faraday / physical-
Faraday residuals are:

| Refinement | Operator | Dual Faraday | Physical Faraday |
|---|---:|---:|---:|
|0|2.72962e-12|1.39337e-10|8.54799e-11|
|1|1.39329e-10|4.95362e-11|3.03465e-11|
|2|6.90051e-11|7.72779e-11|4.73959e-11|

Both residual paths pass simultaneously after balancing their proposed linear
corrections. Full restored Ohm/Faraday and Ampere checks also pass. Source-frozen
replay reproduces every rate and residual trace bitwise, and input coefficients
are unchanged. The independent actual fluid/current/magnetic power contraction
has residual1.68655e-19 in physical power units, or2.69574e-9 of the absolute
sum of its terms. This separate measured cancellation is not a waived solver
bound or an exact-conservation claim.

Independent cylindrical polynomial fields at6/12/24 radial cells reconstruct
analytic B/J/E within2e-10 and satisfy the original rate-work identity within
1e-10. Deliberate iteration exhaustion is atomic. An initial test draft used
a dimensionally wrong E.f work contraction; it was corrected to -Jt.f before
passing, without changing the solver or tolerance. The old successful baseline
is separately checked for bitwise equivalence. Recovery is a linear-iteration
improvement, not a spatial accuracy repair or broad trajectory adoption.
Outputs: `radial16_balanced_probe*`, tests `test_balanced_wall.py`.

## Required next experiments

1. Keep radial16,N128 fixed; next separate axial refinement F32/256,B128/256
   from radial-field refinement F32/128,B256/128, then fluid64. Retain full J/E,
   initial projection, wall/axis, one-sided faces, amplitudes and strict rejection.
   Budget memory explicitly: independent quadrature20 grows with both directions.
2. Localize the remaining large Btheta/vr/Jt rate defect by original physical
   pressure/current/work terms and wall/interface regions. Do not infer a fix
   from conserved total energy or better initial interpolation alone.
3. Once spatial consistency is adequate, repeat full signed100us trajectories
   with genuine time/mesh/residual convergence and original local budgets. Old
   guarded evolution source remainsf255a576; it has not adopted this new solver.
4. Continue the full project: moving heated exact vacuum, transport/ionization,
   coil/end boundaries, prescribed-waveform formation and an experimental UI.
   The user envelope remains~1m length,10–20cm diameter, selectable gas,
   0.1–100mTorr **neutral fill**, and adjustable Helmholtz-like bias. These
   inputs have not been optimized or validated as a laboratory design.

## Completed checkpoint

462 full regression tests pass in230.35s (earlier pre-solver baseline458 in
226.35s). Every calculation, full replay and regression worker is terminal and
collected. Seven verification groups compare every raw array, physical
fingerprint and non-runtime metric, including the explicitly rejected original
fine case. Eleven complete signed rate pairs include the repeated old baseline
and amplitude controls; they are not eleven distinct mesh resolutions.

23 new source archives verify against both their frozen contents and current
source.59 preceding,61 older matched and32 older wall archives remain byte-
unchanged;22 earlier verified source closures still match. Old production
sources and guarded trajectories are unchanged; new solver code is opt-in.
Balanced rate-study source SHA256:
6d1f5b3ff96371e022d486e6acc9ae13c895ae9d155701251a4ad7138b89b7b2.

The complete inventory, test evidence, document/output hashes, unfinished beads
and next action are in outputs/radial16_diagnostics/checkpoint.json. No remote
or publication is configured; the user's modified showcase poster is preserved.
