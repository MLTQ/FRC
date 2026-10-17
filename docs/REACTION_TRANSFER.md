# Reaction return: approximation defect and current candidates

This continues the exact retained-knot rejection described in the earlier
handoff. It is a planar moving-species component study, not a coil-driven FRC
experiment or a recommendation for operating the user's lab. Argon is the
current atomic reference. H2, D2, N2, O2 and other selectable gases require
their own molecular/isotope-specific source and transport physics under
`frc-transport.6`; they are not represented by changing the argon mass.

## Why the trace velocity fails

The actual rejected 20 ns source includes every fine density face and exact
inherited decimal velocity knot. Independent 100-digit elimination of both
the rounded matrix/load and separately accumulated source integrals reproduces
the Ar9+ endpoint magnitude 7,941,858.525 m/s. The rounded solve differs by
1.17e-16 relative; full independent assembly differs by about 3.73e-9 m/s.
The neutral calculation also agrees to about 3.64e-12 m/s. Thus this failure is
an approximation defect despite a well-conditioned diagonally scaled solve.

The physical momentum projection minimizes a density-weighted velocity error.
Very dense material occupying a short part of one continuous velocity element
can determine a slope that extrapolates badly into its dilute remainder. A
separate exact piecewise-constant example has density 1 on [0,epsilon], density
epsilon^6 elsewhere, epsilon=1e-4, velocity 1 only on [epsilon/2,epsilon], and
velocity 0 elsewhere. Its density-weighted P1 endpoints are approximately
-0.25 and 14999.75 m/s although the input velocity lies in [0,1]. The independent
constrained-volume solution is 0.5000375018748124 and -0.2499999943743437 m/s;
the implementation agrees to displayed precision. This static example explains
the mechanism; it is not a physical formation experiment.

Evidence and executable calculations: `outputs/reaction_projection_diagnosis/`.
That directory preserves source closure, exact-input comparisons, raw candidate
states, local errors, independent analytic calculation and inspected plot.

## Candidate contracts and measured limitations

The default remains `projection_method='galerkin'`, `energy_projection='local_total'`.
Two explicit alternatives were implemented without floors, trace cutoffs,
physical-width coarsening, reduced masses or thermal repair.

`bounded_moments` pairs conservative nodal momentum/energy corrections with
convex redistribution. It cures extreme snapshot overshoots but is rejected for
adoption: complete zero-reaction 20 ns controls change final current by 51.25%
and 61.67% at 10/5 ns return intervals. Electric-field changes are 8.99%/7.99%.
This is numerical transfer error despite global energy conservation. The first
warm-decay attempt was interrupted after a nonadvancing clock; subsequent runs
with an explicit clock guard preserve terminal failed prefixes. A charged
attempt also stopped on a JSON serialization error after 10 ns; that reporting
error was corrected before later experiments. These are not completed runs.

`constrained_velocity` minimizes occupied-volume velocity error subject to
exact block momentum and each fine cell's kinetic-energy bound K_c<=E_source,c.
It recovers U_c=E_source,c-K_c once, with no redistribution of cell total energy.
The implemented strictly feasible logarithmic-barrier solve needs a feasible
constant block-mean start; it rejects unsupported states. See its module
companion for equations, tolerances and limitations. Positivity is not accuracy.

On the exact rejected snapshot, Ar9+ endpoints become approximately +/-651.64
m/s (source extrema +/-645.62), with temperature 0.19639–0.28582 eV. Maximum
source-point velocity error is still 197.38 m/s. Neutral minimum temperature is
1.23e-10 eV because a kinetic constraint becomes nearly active. Maximum
normalized momentum error is 6.03e-16 and relative global energy error 3.03e-15;
these conserved totals do not validate the local temperature.

Complete 20 ns zero-source controls preserve the successful retained-knot
comparison: current differences 0.0102096% / 0.0257461%, electric differences
3.14214e-7 / 3.09182e-7 relative for 10/5 ns intervals. There is no local
total-energy redistribution. Direct/10ns/5ns transport counts are 970/971/972.
Evidence: `outputs/reaction_remap_constrained/`; rejected bounded comparison:
`outputs/reaction_remap_bounded/`.

Fixed-input refinement retains every old knot and adds existing density faces
near requested subintervals. Neutral actual velocity cells 38/74/114/138 give
occupied-volume RMS velocity errors 4.8913/4.1505/3.4513/3.0397 m/s. Ar9+ actual
cells 30/58/82/94 give 10.9312/10.5763/8.6317/8.6013 m/s. Ar9+ local normalized
momentum L1 falls 0.12358→0.002496, but maximum neutral temperature change is
0.1050/0.1119/0.1462/0.08875 eV, and Ar9+ maximum source-point velocity error
does not converge monotonically. Therefore locally accurate thermal transfer
is still open. See `outputs/reaction_projection_refinement/`; its first invalid
trace postprocessing attempt is explicitly separated from corrected results.

Additional temperature norms on the same four exact-input levels are saved in
`thermal_refinement.json`. Neutral volume RMS temperature change decreases
0.001956 → 0.001586 → 0.001351 → 0.001027 eV; volume L1 decreases
8.33e-5 → 5.61e-5 → 3.31e-5 → 3.21e-5 eV. Its particle RMS decreases overall
8.68e-5 → 2.48e-5 eV, with a nonmonotone intermediate value. Ar9+ volume RMS
changes 0.0004151 → 0.0004151 → 0.0003225 → 0.0003112 eV. These show improving
average error without established asymptotic order or worst-cell accuracy.
All fine cells and the previous maxima remain included; no small-volume
temperature discrepancy is excluded to improve a metric.

## Trajectories and edge representation

The first constrained, retained-knot 32-cell argon case completes the original
20 ns with 10 ns source intervals: 986 transport steps, 233.27 s elapsed,
21,395 occupied fine cells, minimum width 2.49e-24 m, relative energy change
2.81e-11. Numerical thermal change is -6.10e-12 J/m2; cell total-energy
redistribution is zero. This is useful completion, not convergence or adoption.
The 64-cell case also completes 20 ns: 2,044 transport steps, relative energy
5.60e-12. Its E/J differences from 32 cells are 6.6516% / 32.2137%, and the
electron-front difference is 3.466 mm. These do not establish convergence.
The replacement 5 ns-source experiment retains the full target, as described below.
Evidence: `outputs/reactive_material_argon_constrained_10ns/` and corresponding
`_5ns` / `_n64_10ns` directories. Do not compare unequal accepted times as final
refinement results. The driver's plasma-phase cap also scales with cell count;
that spatial comparison is not an isolated variation of spatial resolution.

The current retained-knot Galerkin warm-decay control itself stalls after three
steps at 74.72 ns of the requested 6.083704 us. New species support endpoints
separated by about 1e-23 m introduce new independent velocity elements and
collapse the acoustic step. This does not invalidate the earlier separately
archived six-case precise-geometry study, but those historical passes are not
passes of today's retained-knot source.

An opt-in edge extension retains every old independent velocity node and
continues its outer affine basis across newly occupied edge material by at most
one adjacent velocity interval. The exact physical endpoint and all fine
density cells remain. Mass, pressure force/work, momentum restriction and
electric force integrate those same shape functions, including negative edge
weights. The convex bounded method cannot use this representation. A separate
source copy allows evaluating it without changing active charged-run sources.
The first complete constrained 16-cell affine case has front error 1.88e-13 m,
velocity error 1.85e-9 m/s, relative energy error 1.28e-13 and population error
8.19e-15. The32-cell case also completes: front error1.12e-14m, relative energy
2.79e-14. The independent16-cell Galerkin control passes. The64-cell constrained
case now completes too, as do all three charged cases: terminal EXIT0,
2701.86s. Finest neutral front error1.42e-15m, relative energy5.68e-14,
population3.33e-14, minimum retained width5.87e-30m. Full raw histories and plot
are preserved and inspected. Bead5.5.7 is closed for this specific defect.

The integrated extension/polynomial/constrained case completes20ns at10ns source
intervals with971steps, relative energy2.84e-11 and no cell-energy redistribution.
E/J change1.53e-6/1.46e-4 relative versus unextended; neutral front changes51.34um
and worst-cell temperatures remain sensitive. The maximum temperature change at
the last source return is0.06161eV versus0.11055eV before extension. These values
do not certify local accuracy. A full5ns extension run is now in progress.

## Equivalent force integration and live refinement

On each existing union interval, E(t)=E0+(E1-E0)t and phi_right(t)=a+b*t for
t in[0,1]. The exact force moment divided by interval width is
a*(E0+E1)/2+b*(E0+2E1)/6; the left moment uses1-a and-b. Precise differences are
formed before float coefficients. This replaces the old degree-two Gauss
integral without changing density cuts, masses or physical approximation order.
Exactly neutral species skip only their identically zero electric-force term.

Full32-cell20ns comparison takes144.21s instead of233.27s, both986 actual steps,
with relative E/J differences1.86e-14/1.70e-8. All3 charged-affine cases pass;
finest relative energy error3.27e-15. The separate unextended warm-neutral
failures are still recorded. The code is integrated; timings are observed runs,
not a hardware benchmark. Source and full results: reactive_material_argon_polynomial_10ns,
precise_material_polynomial, and reactive_material_argon_constrained_comparison.

The old slow5ns run was deliberately interrupted at10ns coupled acceptance;
termination.json records tool exit1/KeyboardInterrupt and the successor. It is
not a physical failure or a complete20ns run. The polynomial and extended 5 ns replays were also retired after diagnosing the
newborn mesh defect below. Both accepted only 10 ns coupled time. Their explicit
termination records preserve that limitation and identify the replacement. No
shortened target is called a success; every active source closure stays frozen.

The specific extreme trace-velocity bug `frc-transport.5.5.4` is closed on the
independent diagnosis, corrected explicit transfer and measured local errors.
Keep `frc-transport.5.5.5/.6/.8`, their parents and `frc-goal` open until their
actual accuracy and full-duration criteria are met. No unit tests were run.

## A distinct newborn velocity-mesh defect

The extended5ns case also slows in its third interval. An exact replay of the
saved10ns endpoint identifies Ar9+ as the limiting acoustic species. Its newly
constructed32-cell nearest-face velocity grid contains a1.2672e-9m interval,
although the nominal cosine grid has minimum spacing4.1755e-5m. Different
nominal targets landed on almost coincident faces within one density cluster.
After100 actual transport steps, this interval contracts to3.3456e-12m and
maximum Ar9+ temperature grows to11.12eV; the actual step is7.51e-15s and only
0.08465ns of the5ns transport interval has advanced. Full split-state and step
records are in reactive_material_argon_extended_5ns/diagnostic_transport_prefix.
This is a diagnostic prefix, not a completed coupled interval.

Bead5.5.8 tracks initial velocity-space selection. An isolated candidate keeps
all inherited independent knots, all density cells and physical endpoints, but
rejects nearly coincident candidate knots before creating a newborn space.
Its explicit spacing fraction is0.5 of associated nominal target spacing;
actual newborn velocity dimension can decrease and is reported. Ar9+ has22
well-spaced cells, minimum width0.26mm, while all21207 density cells remain at
10ns. The full 20 ns, 5 ns-interval run has now accepted 15 ns coupled time,
with 325/203/223 transport steps in the first three intervals and 227,839 density
cells. The last interval is active. This is progress, not yet full-run success or
a local accuracy claim. Source copy:/private/tmp/frc-newborn-mesh-20261016; output:
reactive_material_argon_newborn_5ns. Preserve that imported source while live.


## Optional high-precision arithmetic and integrated checks

`FRC_GEOMETRY_BACKEND=mpfr` selects 192-bit MPFR through the pinned gmpy2 2.3.2
optional dependency. The default is still 50-digit Decimal. Both scalar
construction and local arithmetic contexts set precision explicitly; no global
53-bit default is used. Physical coefficients remain float64 with the existing
species scaling. Every precise coordinate is saved as round-trippable text.
The source identity records arithmetic, precision and dependency versions.

The complete matched 10 ns-source, 20 ns-duration experiment takes 75.43 s
versus 135.89 s with Decimal. All non-coordinate saved arrays are identical,
including populations, velocities, thermal coefficients and topology indices.
Coordinates differ around 1e-50 m. Pointwise temperature maxima across these
slightly shifted piecewise-constant boundaries can remain finite; the largest
volume RMS temperature difference is 2.53e-25 eV. No narrow sliver is deleted
from comparisons. Timing is observational, not a controlled hardware benchmark.

The newborn MPFR 10 ns-source case also completes 20 ns in 971 steps, 72.34 s,
with relative energy error 2.84e-11. Its E/J differ from the old newborn selection
by 6.63e-15 / 6.10e-9 relative. The newly constructed Ar9+ grid has 24 cells in
this trajectory. The integrated main implementation then completes the same
calculation in 74.20 s with identical saved physical arrays and zero E/J
comparison differences. Outputs are `reactive_material_argon_mpfr_10ns`,
`reactive_material_argon_newborn_mpfr_10ns`, and
`reactive_material_argon_integrated_mpfr_10ns`. Main source SHA:
`6148d66b931a92019ea5e14fa7c74a6ea16e945ef16e1dd67f2f72c97f369d3a`.

The newborn spacing option separately completes both full 16-cell analytic
families in `precise_material_newborn`: neutral front error 1.88e-13 m,
population error 8.19e-15 and energy error 1.28e-13 relative; charged front error
2.01e-13 m and energy error 7.60e-13 relative. The inspected figures retain the
independent analytic comparison. The complete MPFR 16/32/64 family and both
full 5 ns-source arithmetic comparisons are still active; these limits must
remain visible until their actual endpoints are collected.
