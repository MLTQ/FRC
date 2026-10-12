# Refined amplitude and residual controls

All workers are terminal and collected, after cd66ec5. Four accepted amplitudes
and two strict rejections have complete fresh protocol replays. Four component,
two spectrum and two magnetic-curl audits also replay completely. No unit tests
were added or run. The full FRC goal and all12 acceptance beads remain open.

## Question and numerical scope

Separate axial B256/320 and radial B320/256 field refinements improved different
current components. This continuation asks how those measurements change with
perturbation amplitude and with the original-equation acceptance tolerance.
Fluid resolution stays64 radial elements /256 axial points. The independently
derived continuum mode stays radial16/N128 with growth12307.15614451037/s.
This is a filled, axially periodic contact fixture with strong radial-wall
constraints, not a formation experiment or a model of physical tube ends.

For each epsilon, evaluate both full nonlinear physical rates and form

    D_epsilon = [F_phys(y0 + epsilon*d) - F_phys(y0 - epsilon*d)]/(2*epsilon).

The physical fluid rate includes the state-dependent quotient derivative when
converting conserved variables; differencing raw conserved coefficients alone
is insufficient. Electric fields are compared with the independent electric
perturbation; density, velocity, species entropy, magnetic and current rates
use the corresponding continuum time derivatives. For every scalar component
and full vector, errors and actual changes are normalized by that component's
independent continuum L2 norm with volume measure2*pi*r*dr*dz.

The three amplitudes1e-5,1e-4,1e-3 separate sensitivity from mesh error. A
centered derivative of a smooth exact rate has an O(epsilon^2) truncation term;
finite solver and arithmetic errors can be amplified by division by epsilon.
These observations do not identify a particular error source or prove that
the middle amplitude is an asymptotically valid linear-response limit.

Every accepted rate must satisfy the unchanged1e-10 original residual bound.
Separate1e-11 controls use epsilon1e-4 and all existing iteration/refinement
budgets. A strict failure retains its exact saved inputs and every failed
attempt; it has no accepted physical rate. Its full fresh replay verifies the
rejection protocol, not physical accuracy. Resource interruptions are distinct.

The original component-stress physical source remains
000fd344b524568a655485308fb72d754aa8c144d26cf3a8e04c4d947ffb53ef.
Only explicit orchestration/reporting modules are new. No equation, boundary,
time integrator, diffusion, filter, atmosphere, mass, floor, state or energy
repair is introduced. There is no time evolution in these signed-rate controls.

## Complete physical runs

All four new amplitude primaries and fresh full replays pass bitwise. At q20,
the existing whole-vector current-rate errors are:

| Field radial / axial | epsilon1e-5 | epsilon1e-4 baseline | epsilon1e-3 |
|---|---:|---:|---:|
|256 /320|0.0811128%|0.0801366%|0.0810908%|
|320 /256|0.0563520%|0.0552027%|0.0566025%|

These norms hide weaker current components; the independent component and
actual magnetic-rate curl comparisons below expose them. Worst residuals for the four
new amplitudes are7.34506e-11,2.01063e-11,9.95202e-11,2.28167e-11. The radial
small-amplitude case is close to, but below, the unchanged1e-10 bound.

Axial1e-11 strictly rejects on the negative sign after all five outer attempts
(inner tolerances1e-11 through1e-15). At the last attempt its operator,
Faraday-dual and physical-Faraday residuals are3.57503e-11,6.14598e-11 and
1.65024e-11. The radial1e-11 control also rejects on the negative sign after
all five outer attempts; its last operator, Faraday-dual and physical-Faraday
residuals are4.70168e-11,8.29280e-11 and1.92625e-11. Both entire rejection
protocols replay bitwise, including their saved inputs and complete messages.
All five saved arrays in each rejected case have identical shapes, dtypes and
bytes to the same-grid accepted baseline. This establishes saved-input equality,
not post-call in-memory immutability during the rejected derivative.
No partial iterate is used for physical norms. These runs reproduce the strict
failure on both refined grids; they do not independently identify its cause.

Every physical worker is terminal and collected. Accepted amplitude
primary/replay peak RSS in GiB is:

| Grid / amplitude | Primary | Complete replay |
|---|---:|---:|
|256 /320,1e-5|7.174|6.681|
|256 /320,1e-3|6.830|6.787|
|320 /256,1e-5|6.809|7.574|
|320 /256,1e-3|7.489|7.437|

Axial tighter-bound primary/replay peak6.672/6.549GiB; radial6.798/6.898GiB.
No resource guard fired, and final extrema also satisfy8GiB/256MiB. The least
free disk observed by any physical worker is605.469MiB. Observed peaks are not
guaranteed memory bounds for a further or jointly refined grid. Verified
byte-identical completed replay files share immutable storage with all paths
and hashes retained in the campaign manifest.

For the accepted new amplitudes, instantaneous
abs(P_fluid+P_current+P_magnetic)/sum(abs(P)) is at most3.98223e-8. This is a
relative power-exchange diagnostic at two states, not local charge/enthalpy/
Poynting closure, a time-integrated conservation bound or a formation result.

## Independent physical component errors

Every entry below is a q20 physical L2 rate error in percent, normalized by
the component's own continuum norm. All component audits reproduce the original
physical/reference/weight fingerprints and replay every saved array and metric.

| Field radial / axial | epsilon | Jr | Jtheta | Jz |
|---|---:|---:|---:|---:|
|256 /320|1e-5|4.17506|0.0157584|1.10491|
|256 /320|1e-4|1.76199|0.0157584|1.09567|
|256 /320|1e-3|1.77354|0.0157584|1.10924|
|320 /256|1e-5|6.70354|0.0121202|0.746810|
|320 /256|1e-4|5.90438|0.0121202|0.734849|
|320 /256|1e-3|5.90970|0.0121201|0.755276|

The actual small/baseline radial-current differences are3.83376% and3.19588%,
although the full-current vector norms change little. Of their squared norms,
90.2292% and91.5176% lie in the interior5%-to95%-radius band; only9.59082%
and8.28437% lie in the outer5%. The sensitivity is not confined to the wall.
These bands mask original Gauss nodes and do not split cells at their edges.

Larger amplitude is not an accuracy repair: it slightly increases radial/axial
current and total-current errors, and also increases radial-velocity and
azimuthal-field errors. Azimuthal-current errors barely change. The intermediate
amplitude performs best for those former channels among these three samples;
this is not an amplitude-convergence proof or a selected laboratory waveform.

## Cell spectra and actual magnetic response

All120 direct squared integrals per grid agree bitwise with the separate
component audit,240 in total. The radial-current error above cell mode|k|64,
normalized by the same continuum component, is:

| Grid | epsilon1e-4 RMS | epsilon1e-5 RMS | Small / baseline | High-band share of squared difference |
|---|---:|---:|---:|---:|
|256 /320|0.00368476|0.0368299|9.99520|93.7982%|
|320 /256|0.00305835|0.0298545|9.76164|87.2576%|

These RMS entries are fractions, not percentages. Multiplication by epsilon
gives nearly constant high-band scales on each grid. This is consistent with
approximately amplitude-independent numerical noise amplified by differencing,
but does not identify its source. The radial-grid baseline still has dominant
lower-band error. Cell-translation spectra retain intra-cell coordinates and
fold aliases separated by the axial cell count; no mode or state is filtered.

Independently reconstruct S=curl(Bt)/mu0 inside the polynomial cells:

| Grid / epsilon | Jr error | Sr error | Jr-Sr norm |
|---|---:|---:|---:|
|256 /320,1e-5|4.17506%|3.19274%|1.32585%|
|256 /320,1e-4|1.76199%|1.73644%|0.131541%|
|256 /320,1e-3|1.77354%|1.76837%|0.0151631%|
|320 /256,1e-5|6.70354%|6.27542%|1.14822%|
|320 /256,1e-4|5.90438%|5.87664%|0.119842%|
|320 /256,1e-3|5.90970%|5.88662%|0.0297513%|

Thus the sensitivity also occurs in the differentiated magnetic response;
current reconstruction alone is not a demonstrated cure. Strong curl omits
distributional face jumps, so J-S is not itself a failed weak Ampere equation.
All48 current error/reference integrals across the four controls and two
quadratures match the separate component audits bitwise.

## Quadrature, resources and preservation

Whole-volume component relative norms agree between q16/q20 within3.43905e-15
absolute; original full-rate norms within4.29778e-16; current/strong-curl norms
within2.68368e-16. Band-norm differences are larger: axis1.58380e-5,
wall4.96785e-5 and interior1.31929e-6. Do not substitute volume agreement for
regional accuracy. The largest spectral-fraction difference is2.14728e-6,
for the precision-scale density amplitude difference. Parseval defects are
at most3.14934e-16; curl squared-error decompositions close within2.29706e-16.

The independent-audit process peaks at1.464GiB; its minimum observed disk
headroom is294.531MiB. Those are stage-boundary observations/preflights, not
continuous guards. Seven further verified replay NPZs share immutable storage
after the earlier one-file sharing, preserving every path, byte and hash.
No interruption occurred in this continuation.

Scientific PNG/PDF: outputs/refined_control_report. The PNG was visually
inspected; the PDF uses the same figure but was not separately inspected.
The report retains complete component/region/curl/spectral metrics, original
residuals, work diagnostics, resource observations and input/artifact hashes.

Four new modules have companions. All58 new source archives,446 preceding
archives,164 preceding sources,228 preceding checkpoint evidence files and22
active older source closures verify. The user-modified poster is preserved.
Final168-file source union:
5deb5dc6c4a3229d60d13d33b2a182716aadd61fcd1fcfa1fa5264ef6dff4966.
Final inventory: outputs/refined_control_diagnostics/checkpoint.json.
Four accepted simulation pairs, two rejection-protocol pairs and eight audit
pairs make14 complete replay groups and14 raw-file comparisons.

## Next concrete work

Resume the pending radial1/N64 full signed100us trajectory comparison on
F8/64,B32/64 against preserved F4/64,B16/64, with60/120-step controls, original
stage/endpoint/induction tolerances, independent field/conservation audits and
full stage/state replays. It is a different reference from this radial16 study;
its near-grid growing response and time convergence remain unresolved. Budget
and guard actual memory/disk; source-derived retained modal storage is only a
lower estimate. Read JOINT_E2E_DESIGN and the split_mode_evolution companion.

Keep refined-grid amplitude noise and1e-11 rejection open. Any arithmetic or
expanded electric representation must preserve the real equations, force/work
adjoints and original bounds, with full physical controls before adoption.
Joint B320/320 refinement also remains unbudgeted; do not infer feasibility or
accuracy from independent direction increments. Local budgets, heated moving
exact vacuum, transport/ionization/validity, real prescribed coils/ends,
converged formation and UI remain required. No3-D stability/confinement claim.
