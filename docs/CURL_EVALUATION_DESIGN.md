# Projected-curl arithmetic and physical cell-mode errors

Continuation from398f9a3, classified as verified progress. Keep the complete
experimental sandbox goal open. No unit tests are added or run.

## Captured rejected state

For the unchanged raw curl C and wall projection P, compare F e with P(C e),
where the stored F is the binary64 matrix product P C. Preserve the originally
evaluated residual vectors at every captured iterate. Rebuild the original
dual lift D = H T and require bitwise RHS, rate, dual-Faraday and electric
distributivity vectors before interpreting changes in operation order.

Measure the assembled, sparse-factored, dense-radial-factored and paired-stencil
actions on e0, ei, rounded e=e0+ei and its independently verified addition low
part. At the final iterate use exact Decimal products and short raw CSR sums,
rounded once to binary64, as an independent raw-curl comparison. Projection
and D remain binary64. No alternative action accepts a rate or replaces a
solver check. A smaller diagnostic residual does not demonstrate a repair.

## Physical component errors

Reconstruct full signed responses at the same original quadrature points,
verifying original physical fingerprints. Fourier transform errors only over
axial-cell index at each radial/intra-cell point. Discrete Parseval provides
an independent conservation check on the spectral partition. Retain each
component's own full reference norm, complete signed-mode spectrum, quadrature
differences and full original/replay evidence. These cell modes fold continuous
harmonics differing by the number of axial cells; do not claim an unaliased
continuous spectrum. No accepted field is filtered.

The complementary current-rate audit compares Jt with the independently
differentiated strong curl of the actual centered Bt, inside every polynomial
cell. With S=curl(Bt)/mu0 and independent target Jref, record the exact vector
identity Jt-Jref = (S-Jref)+(Jt-S), its squared norms and signed cross term.
This distinguishes magnetic-response error from weak/strong reconstruction
differences without calling Jt-S a violated weak Ampere equation. Distributional
face jumps are excluded from this volume comparison and remain separate gates.

## Additional full amplitude control

The original complete spectra show that amplitude sensitivity is mostly above
the independent reference's axial band, whereas the larger-amplitude error is
mostly inside it. Add epsilon1e-3 to the unchanged F64/256,B256/256 full signed
workflow, alongside the preserved epsilon1e-4/1e-5 cases. Keep tolerance1e-10
and every physical term, original-equation check and independent q16/20 audit.
This is a centered-derivative control, not a proposed physical experiment or
permission to increase amplitudes until a bound passes. Independently replay
the complete positive/negative computation, then reconstruct individual current
errors and cell spectra. A diminishing high-mode contribution with unchanged
low-mode error would distinguish amplitude sensitivity from spatial error;
nonlinear truncation must remain visible and is not assumed negligible.

## Unchanged physical scope

Original filled-periodic argon fixture, physical masses, actual density and
impermeable wall remain unchanged. No numerical diffusion, atmosphere, floors,
physical transport or energy repair is added. Actual current-component
convergence, full signed trajectories/local budgets, heated moving exact vacuum,
transport/ionization, real coils/ends, formation and experimental UI remain open.
Large geometry builds/replays run in separate processes with resource guards.
