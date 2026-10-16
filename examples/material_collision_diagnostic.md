# Collision sensitivity of an actual completed pulse

Evaluate the finite-drift paired Maxwellian collision moments at each recorded
state of a completed configured experiment. Use the source's declared fixed
logarithm for a paired run, or explicit exploratory log10 for prescribed eta;
the screening/cutoff prescription is not established for either. Compare eta*J²
where active, small-drift Maxwellian drag using the same
logarithm, finite-drift drag and separate electron/ion heat. No field, temperature
or current is changed, and no candidate trajectory is evolved.

The integral along this **original** trajectory is a sensitivity diagnostic;
it is not predicted heating of a future pulse using altered transport. Its sparse
summary-state quadrature is also distinct from the original solver's every-step
Joule integral. Each sample's heat sum and drag are retained without correction.
Maximum local relative-drift relaxation rate times original dt screens one
source timescale; it is not a full-system stability criterion or timestep choice.

The example freezes its own source/companions and hashes input states/summary,
then emits two PNG/PDF figures plus actual values. The model-validity limitations in
`maxwellian_pair_collisions.md` remain. This is physical model exploration, not
a pass/fail suite, convergence certificate or new verification framework.

`regime_content` extends the existing NRL thermal diagnostics to every actual
positive plasma quadrature point. It records fractions of particles, volume,
electron thermal energy and candidate drag power in explicitly overlapping
review regions. Unavailable thermal coefficients remain unavailable; no log
floor supplies a replacement. The hypothetical path-energy fractions use the
same unchanged-state time quadrature as the original collision sensitivity.

`density_log_gradient` differentiates the configured profile analytically inside
each actual Q1 material cell. For the meridional map A=dx/da, R its physical
radius and initial label r, det(F)=(R/r)det(A). Thus

    grad_x ln(n) = A^-T [grad_a ln(n0) - grad_a ln(det(F))].

The known n0=nc*[(1-r²/R0²)(1-z²/L0²)]^(3/2) supplies its exact reference
gradient. The bilinear map's mixed derivative supplies grad(det(A)); rotation
does not change this axisymmetric scalar gradient or determinant. No cloud
differences, smoothing or density floor enter the diagnostic. Derivative jumps
across cells remain; temperature and field gradients are not implied by Ln.

Thermal mean-free-path ratios use Ln=1/|grad ln(n)| and its direction along actual
B. Electron weak/strong-field groups are separated by omega_ce/nu_ei=1. Ratios
>=0.1 and logs<5 are existing conservative review thresholds, not physical
transitions or a proof of invalidity/validity. Ion isotropic density ratios are
reported without asserting a complete magnetized transport criterion. The
electron gyro/Debye review is retained in JSON. Reference: NRL Formulary (2009),
printed pp.34,37–39; `docs/TRANSPORT_REGIME.md` retains the source link and formulas.
Large drift and non-Maxwellian distributions limit these thermal estimates.
Two scientific figures separate collision-power sensitivity from physical-region
content; their values are saved before rendering. No trajectory is modified.

For a completed paired-Maxwellian source run, use its declared logarithm and
label the reconstructed pair power as the actual implemented collision law.
For a prescribed-eta source, retain hypothetical log10 sensitivity. The role is
explicit in JSON and figure labels; legacy `candidate_*` field names remain for
compatibility. The unused zero prescribed-eta curve is omitted for paired runs.
The original every-step total heat is retained beside the sparse path integrals.
An actual implemented law still does not establish kinetic or transport validity.
