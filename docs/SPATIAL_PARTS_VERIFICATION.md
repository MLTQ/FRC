# Spatial contributions and local current flux — 2026-10-12

This continuation made verified progress from fa597fc. All physical workers
are terminal and their actual exits were collected. No unit tests were added
or run. The full FRC goal and all12 acceptance beads remain open.
The original physical solver, states, boundaries and tolerances are unchanged.

## Complete original-equation decomposition

All12 signed states from the tangent probe reproduce the original complete
fluid, B_t, J_t and E arrays bitwise. Ten groups retain bulk, ideal/Lorentz,
Hall, changing-density inertia, advection, stretching/drift stress, self drift,
electron/ion pressure-entropy and wall reaction. The bulk wall reaction is
also retained separately before the declared inactive-slot constraint.

All96 separated induction solves pass the original1e-10 bound. Maximum
component residual is1.021886e-11; species-potential-restored Ohm backward
error is at most3.079e-15. Recombined Ohm/Ampere/Faraday backward errors are
at most2.360e-16/8.759e-13/7.811e-15. Complete native fluid sums are bitwise
identical to the full solver. No term was removed or used to evolve a state.

Independent solves do not produce an exactly linear numerical map: maximum
single-sign native B/J sum differences are1.18844e-5/2.56069e-4 relative.
These remain exposed despite passing backward errors. In the independent
physical centered responses, maximum component sum differences are1.33244e-7
for the saved error and3.95787e-5 for the known mode, both at h=1e-6.
The latter is axial-current response; small-signal subtraction matters.

All12 raw signed contribution records and non-runtime metrics replay exactly.
The physical audit reproduces every original full-field fingerprint at q16/20
and itself replays all metrics and its raw integral table bitwise. Source and
raw evidence: `trajectory_spatial_parts`, its `_replay_verification`,
`trajectory_spatial_parts_audit` and its `_replay_verification` under outputs.

## Physical directional contributions

At h=1e-5 and q20, contributions to the Btheta directional projection are:

| Term | Saved error (/s) | Independent mode (/s) |
|---|---:|---:|
| Ideal / Lorentz | -534.112 | -252852.851 |
| Hall | -664586.832 | -652684.554 |
| Electron pressure/entropy | +1615416.282 | +914289.969 |
| Ion pressure/entropy | +22.4433 | +12.3500 |
| Complete sum, including other terms | +950317.782 | +8764.91350 |

These are signed projections in each physical component's own norm, not
eigenvalues, energy fractions or probabilities of a cause. For saved-error
electron entropy, bulk contributes-44.9340/s and electron current transport
950243.810/s. Radial velocity has-93770.468/s bulk pressure,1045179.006/s
Lorentz and328.551/s wall reaction. The coordinate metric excludes that
inactive wall slot; the physical reconstruction retains the reaction.

The changing-density/current-advection/stretching/self-drift projections shrink
as h². Their adjacent-amplitude relative differences approach0.99 because
their linear limits vanish, not because they fail to converge. Their actual
norms and full-response-relative sizes remain in the audit; they are not
discarded. Significant signed contributions retain separate amplitude records.

## Independent continuum and interior comparisons

A separate Fourier/Bessel reference splits original ideal, Hall and species
induction sources with physical masses and analytic pressure curl. Its summed
rates reproduce the frozen continuum mode within1.126e-8 relative in the worst
component (radial current); all fluid components are within3.239e-11. These are
reference consistency checks at N64, not a new continuum refinement study.

Known-mode h1e-5 axial-current errors are4.56401e-4 ideal,1.68286e-4 Hall and
2.67214e-4 electron, each relative to that term's own continuum norm. The full
current response error remains5.08482e-3 against its much smaller net target.
Electron entropy's separate continuum error is4.65071e-6. Large opposing
responses make the net response more sensitive to component errors.

For the actual saved error, bulk density and electron-entropy discrepancies
against the independently reconstructed strong interior equations are0.552022
and0.550012. Bulk radial-acceleration discrepancy is1.010913. These omit
distributional interface and boundary terms; they cannot alone establish a
discrete inconsistency. Adding the volume-only conservative or half-divergence
entropy corrections does not restore agreement: electron-entropy discrepancies
become1.065565/1.095158. Neither alternative is inserted into the solver.

Across all stored q16/20 integral entries, maximum squared-norm relative change
is1.29021e-10 (a vanishing nonlinear axial-velocity term with integral1.72e-20).
Maximum signed-product relative change is1.94517e-8 (a tiny wall-reaction
product). Both full tables and the locations/values of these extrema are saved
in `trajectory_spatial_parts_diagnostics/quadrature_comparison.json`.

## Independent one-sided physical charge audit

The charge audit uses explicit-cell root-product edge bases, full physical
volume derivatives and both normal traces at every field-cell face. It also
checks analytic Jr=r*k*cos(k*z), Jz=-2*sin(k*z), whose divergence is exactly zero.
All4 q16/20 direction records, trace arrays, divergence fields, cell integrals
and metrics replay bitwise. The previous physical divergence integrals match
exactly; this is independent trace evaluation, not a different current.

| Dimensionless measure, q20 | Saved error | Independent mode |
|---|---:|---:|
| R ||div J|| / ||J|| | 0.5024685445 | 0.000377938468 |
| sqrt(R) ||radial normal jump|| / ||J|| | 0.0560483418 | 2.06491e-8 |
| sqrt(R) ||axial normal jump|| / ||J|| | 7.33080e-5 | 0.000118298134 |
| sqrt(R) ||wall normal current|| / ||J|| | 2.88224e-17 | 2.28927e-16 |

Volume and surface norms have different dimensions before this explicit
normalization. Direction amplitudes are arbitrary fixed normalization factors,
not physical trajectory current amplitudes. Both quadratures agree closely.

Every cell's integrated physical divergence agrees with its one-sided outward
current flux to at most5.38902e-16 relative to absolute face fluxes; maximum
cell defect is4.59910e-14A for the normalized saved-error direction. The analytic
zero-divergence control's relative cell-flux defect is at most1.057e-16.
Global volume divergence and net boundary flux are near cancellation roundoff,
while local divergence and internal radial jumps are substantial. This verifies
the measurement and locates an unresolved local charge-consistency problem;
it does not equate the solver's conserved graph charge with geometric fluxes.

## Resources, presentation and next action

All four guarded contribution/physical workflows pass8GiB/256MiB including final
observations; maximum peak4087185408 bytes, minimum free6225920000 bytes.
Charge audit/replay stage observations give maximum429506560 bytes and minimum
free6144000000 bytes; those disk observations are not continuous monitoring.
Every actual worker exit is0, including both complete replay verifications.

Inspected report: `outputs/spatial_parts_figure_v2/verification.png`; PDF is the
same figure and was not separately inspected. The first two presentation
variants and sources remain preserved, including their observed layout issues.

NEXT derive and independently measure the complete weak species entropy flux
through these internal interfaces, retaining physical current jumps, potential
projection/adjoint, actual fluid mass metric and wall reactions. Check whether
the fast direction occupies poorly represented potential/current subspaces;
compare a field-refined control only with explicit field/current transfer
errors. Do not remove electron/Hall terms or add filtering as a fix. Causality,
full signed/time/arithmetic trajectories, local energy/enthalpy/Poynting budgets,
heated exact vacuum, transport/ionization/validity, real coil/end coupling,
converged formation and the experimental UI remain unfinished.
