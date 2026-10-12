# Curl evaluation and physical current-component evidence

The original experimental sandbox goal remains open. No unit tests were added
or run. These are filled-periodic verification fixtures, not selected lab
conditions. Equations, numerical scope and comparison protocol are in
CURL_EVALUATION_DESIGN.md. All workers are terminal and collected.
Five complete verification groups
match every saved array and non-runtime metric. These include four physical
or arithmetic diagnostics and one full signed instantaneous rate calculation;
they are not five formation experiments or full time trajectories.

## Original rejected curl actions

All four captured original RHS, magnetic rate, dual-Faraday and electric-curl
distribution vectors reproduce bitwise from the rebuilt original geometry and
dual lift. The complete workflow and its fresh standalone replay match every
saved array and non-runtime metric. Original/replay elapsed times are 223.731 /
211.941s. No solver or accepted state was changed.

For the final rejected iterate, norms divided by the original RHS norm are:

| Evaluation | Electric distributivity | After addition-low accounting | Dual Faraday |
|---|---:|---:|---:|
| Assembled wall matrix | 6.87232e-11 | 4.55034e-11 | 5.88820e-11 |
| Curl then projection | 6.75088e-11 | 4.49539e-11 | 6.45170e-11 |
| Subtraction-first raw curl | 5.02775e-11 | 1.25835e-12 | 5.61985e-11 |
| Exact raw products/sums, then projection | 5.02657e-11 | 2.96583e-15 | 5.61939e-11 |

Simple factoring does not fix the discrepancy. Independent Decimal arithmetic
computes each raw CSR row from exact binary64 operands at 2200-digit precision,
then rounds once. Its scope is the full final-iterate e0, ei, e and addition-low
arrays. Projection remains binary64. The original raw CSR full-electric error
against that reference is 3.16250e-11/RHS; subtraction-first reduces it to
8.99408e-13/RHS, after the same projection and dual lift.

Accounting for lost electric-addition bits plus exact raw arithmetic nearly
restores distributivity, but the actual saved iterate still violates 1e-11.
This comparison neither accepts a failed derivative nor verifies a changed
solver. Any future arithmetic implementation must preserve the original real
equations, adjoints/work, physical fields, all bounds and full signed replays.

## Physical axial cell-mode partition

The complete epsilon1e-4/1e-5 component reconstruction and its replay match all
original field/reference/weight fingerprints. All 120 component squared
integrals are bitwise equal to the preceding physical component audit.
Cell-Fourier Parseval defects are at most 2.03e-16 at quadrature20.

For radial current, 86.672% of squared amplitude difference lies above |k|64;
52.260% is in |k|>=96. By contrast, only 0.2195% of squared error at epsilon1e-4
lies above |k|64. Its largest error mode is |k|20, while the reference is
dominated by |k|2. At epsilon1e-5 the high-band error fraction rises to 17.448%.
Thus the amplitude sensitivity and dominant larger-amplitude spatial error
have different spectra. These are axial-cell translation modes with the full
intra-cell coordinate retained; continuous harmonics separated by 256 fold
together. No accepted state or field was filtered.

## Independent magnetic/current-rate comparison

The complete two-amplitude audit and fresh replay agree bitwise. Independently
differentiate the actual centered magnetic rate inside polynomial cells to
form S=curl(Bt)/mu0. Compare J, S and the independent current-rate reference:

| Amplitude / component | J error | S error | J-S norm |
|---|---:|---:|---:|
| 1e-4 radial | 0.05969389 | 0.05942468 | 0.001107695 |
| 1e-5 radial | 0.06672077 | 0.06293120 | 0.01077299 |
| 1e-4 axial | 0.01187001 | 0.01003199 | 0.003452958 |
| 1e-5 axial | 0.01193623 | 0.01009534 | 0.003461323 |

All are physical volume L2 norms divided by each component's own reference.
The signed cross term closes the squared-error decomposition to 2.16e-16
relative. All original current/reference/weight fingerprints reproduce. The
radial discrepancy remains in the differentiated magnetic response; changing
the current reconstruction alone is therefore not a demonstrated cure.
Strong curl excludes distributional face jumps. J-S is not itself a violated
weak Ampere equation or proof of a defective discrete adjoint.

## Additional full amplitude control

The unchanged F64/256,B256/256 full signed epsilon1e-3 workflow completed in
477.099s. Both signs pass every original check; worst residual 1.91046e-11 is
below the unchanged 1e-10 bound. Independent q16/20 audits give radial-velocity,
toroidal-B and total-current errors 0.001893220  / 0.005643362  / 0.000887386.
These exceed the corresponding epsilon1e-4 values. The larger signal is an
amplitude control, not an accuracy repair or a proposed experiment setting.
Complete physical component-spectrum replay and complete signed rate replay
are bitwise verified. The full rate replay took 488.773s. Individual
radial-current errors are:

| Amplitude | Total error | Within reference band | Above reference band |
|---|---:|---:|---:|
| 1e-3 | 0.05975116 | 0.05975053 | 0.0002750160 |
| 1e-4 | 0.05969389 | 0.05962835 | 0.002796566 |
| 1e-5 | 0.06672077 | 0.06062118 | 0.02786995 |

Band norms combine in quadrature. Multiplying the high-band norm by amplitude
gives 2.75016e-7, 2.79657e-7, 2.78699e-7. This near-constant scale is consistent
with amplitude-independent numerical noise amplified by the centered
derivative; it does not identify the responsible arithmetic operation. The
roughly 6% in-band error persists. Larger amplitude changes velocity/toroidal-B
responses as well, so it is not a substitute for spatial convergence.

Across all component spectra, q16/20 physical relative norms differ by at most
3.44e-15 absolute. The largest mode-fraction change is 3.38e-6, for a density
amplitude difference already near floating precision. The independent
magnetic/current-curl norm changes are at most 1.77e-16; the new full signed
workflow's volume-error changes are at most 2.84e-16. Complete values and
scales are retained in outputs/curl_evaluation_diagnostics/quadrature_metrics.json.

Scientific figure: outputs/curl_spectrum_report_layout/verification.png,
visually inspected. The earlier layout and its frozen source are also
preserved; the layout-only wrapper changes no values or physics.

## Remaining work

Resolve the dominant spatial magnetic/current-component error with directional
refinement and independent component norms. An accurate curl-action/expanded
electric representation can be investigated, but no current result proves a
full strict-tolerance repair. Then complete signed trajectory convergence and
local budgets, moving exact vacuum, transport/ionization, real prescribed
coils/ends, converged formation and the experimental UI. Full 3-D stability and
confinement remain separate claims, with no evidence supplied by these fixtures.

## Preservation and continuation

All six new modules have companions and are frozen. Eighteen new source
archives, 300 preceding archives and 22 earlier source closures are verified.
Final source union: 140 files, SHA
b0ee5f17e5baf0fa4fc564d3fe40437173b05479d16210e2ef040f0a1d0ce922.
Complete input/source/replay/worker inventory:
outputs/curl_evaluation_diagnostics/checkpoint.json. No resource guard fired,
no prior calculation source changed, and no unit tests were added or run.
The user's poster is preserved outside the research commit. Research and
beads remain local; no remote or publication. All 12 acceptance beads and the
overarching FRC goal remain open.
