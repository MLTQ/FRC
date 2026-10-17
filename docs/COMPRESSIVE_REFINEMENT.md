# Compression refinement: completed physical comparisons

These are planar component experiments, separate from formation and the reacting
precise-coordinate solver. Physical Ar-ion/electron masses, all material and
exact empty exterior are retained. No unit suite was run. Neither conservation
nor agreement with a numerical reference establishes lab accuracy.

## Method and boundaries

`plasma_bench/compressive_refinement.py` selects unresolved compressive velocity
jumps using the existing minmod residual relative to actual sound speed. The
explicit threshold and minimum child width are numerical resolution choices.
Bisect selected full cells, retain every old node and interpolate P1 velocity.
Constant children halve N/U, preserving represented P0 density/thermal fields.
The separate limited-linear option reconstructs bounded pressure before halving:
Uleft/right=U/2-/+s*h²/[8(gamma-1)], with minmod neighboring center secants s.
It preserves parent total U but deliberately redistributes thermal energy.
Per-event particle/energy/momentum changes and redistribution are recorded.

Compatible Euler-Poisson transport uses gamma5/3 and minmod numerical viscosity
(.5,2), whose mechanical work heats the same material. This is numerical shock
stress, not physical viscosity. Supports move freely inside remote field ends;
no atmosphere, ghost material, remapping, thermal repair or species cutoff.
Capacity or actual thermal/geometry failure terminates the original request.
See helper and driver companions for the equations and saved-state contracts.

## Independent neutral shock

Physical neutral argon occupies±50mm, with left/right n=1e19/1.25e18m^-3 and
T=1/.8eV. Remote ends are±.4m. All cases complete3.860718379770096us; compare
central |x|<25mm to the independent ideal Euler Riemann solution before exterior
rarefactions arrive. The exact shock/contact/fan cuts enter the quadrature.
Threshold .02; minimum child widths50/25um for64/128 initial cells.

| Children | Initial/final cells | Steps | Density L1 | Velocity L1 | Pressure L1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Constant |64/371|1155|.00671549|.0135080|.00972332|
| Constant |128/575|2300|.00397639|.00775253|.00543964|
| Limited linear |64/364|1155|.00673343|.01346865|.00970214|
| Limited linear |128/564|2301|.00385802|.00766029|.00524218|

Norms scale density by1e19m^-3, speed by sqrt(eV/mAr), pressure by1e19eV/m^3,
and integrate over the50mm central interval. Energy drift is4.18588e-7/2.09294e-7
for the64/128 families. Event energy changes are at most6.30e-16 relative.
Limited-linear cumulative absolute thermal redistribution is5.18531e-4/
3.61432e-4J/m²; largest child temperature change is.0164568eV. These measure
numerical redistribution, not added physical heat. Errors decrease, but the
added cells do not demonstrate greater efficiency than fixed finer grids.

The first64-cell constant run completed dynamics and saved endpoint/budgets but
failed JSON serialization of a NumPy boolean. Preserve its reporting_failure.json
and use the fresh full replay for complete event records; no dynamics changed.

## Full physical-mass plasma: both adaptive candidates rejected

Both species initially occupy±10mm at n=1e14m^-3, Ti/Te=.2/2eV, zero velocity
and electric field, cosine labels, remote ends±.2m. Complete the original20ns.
Only electrons refine, threshold .01, minimum child width5um, plasma-phase cap
.01 plus actual acoustic restriction. No reactions/collisions in this component.

| Method | Final electron cells | Steps | Maximum energy drift | J difference | E difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fixed128 |128|existing|historical|34.2024%|2.44725%|
| Constant split |3664|28786|6.96381e-8|34.1688%|.973794%|
| Limited-linear split |2904|26348|7.47130e-8|32.7731%|.911228%|
| Fixed256 |256|existing|historical|19.0707%|1.02490%|
| Fixed512 |512|4860|5.34015e-10|8.07805%|.330304%|

J/E are relative L2 differences from the same full1024-cell numerical endpoint.
Two-point quadrature on the exact floating-coordinate union retains all support
and empty exterior. Time caps differ: fixed grids use .04*32/N. Thus this is
not isolated spatial order. The1024 reference itself remains unconverged.
Constant/linear current plots show new large oscillations around5–11mm in
absolute position. Better E and conserved totals do not make these acceptable
current refinements; neither option is adopted in the reacting/formation solver.

The limited-linear plasma run has1316 refinement events, maximum per-event
energy change4.01e-16, cumulative absolute redistribution1.01932e-8J/m² and
maximum child temperature change.0787898eV. The constant case has1640 events.
All full initial/final states, actual dt/budgets, events and frozen sources are
saved. The driver intentionally avoids padded variable-mesh histories; source,
initial data, actual steps and deterministic refinement rules enable replay.

## Evidence and next action

Constant source SHA065369f976b78d094c6a6d87aee0e2f2e9594d9c164c6fbf917077388e819043.
Limited-linear source SHA6c61fe5394eed753f26d2a2bf0e21644cb022ad4a830c7cea5ab6a4832c114fd.
Output directories start `outputs/compressive_refinement_`; comparison source,
JSON and inspected PNG/PDF are in `outputs/compressive_refinement_current`.
Use `.venv/bin/python -m examples.compressive_refinement` with a fresh directory;
see its companion for exact parameters and selectors.

The full fixed 2,048-cell trajectory completed 20 ns in 19,402 steps (115.005 s),
with all refinement disabled and phase cap .000625. Maximum relative energy
drift is 7.49597e-11, number drift zero, momentum drift 3.46230e-27 kg/(m s).
Source: 1d9e31d9385bb77492bdb08c7fc1b28fe112ab257730c9b061af713432f816d4.
Its small endpoint/budget archive avoids another multi-gigabyte padded history.

Successive fixed-grid current L2 differences are 8.07805% (512 to 1,024) and
5.38767% (1,024 to 2,048); E differences are .330304% and .156849%.
Of the latter squared current difference, 82.3211% lies inside |x|<9 mm.
The inspected fixed-grid plot shows continuing changes in narrow current peaks.
This is improving numerical agreement, not established current convergence;
do not treat the 2,048-cell result as exact or adopt either adaptive method.

Next collect the preserved all-heavy captured 5 ns stage, then evaluate the
full all-heavy 5 ns-source coupled trajectory and source/space sensitivity.
Do not restart or shorten the live stage. The current branch requires a different
transport improvement; further threshold tuning of the rejected split methods
is not the default next action. Preserve the new fixed reference for comparison.
Thermal/current convergence, physical gas closure, boundaries and formation
validation remain required by the overarching goal.
