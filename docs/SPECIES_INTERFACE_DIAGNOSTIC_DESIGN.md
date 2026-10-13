# Next complete weak species-entropy interface diagnostic

Start from the completed SPATIAL_PARTS_VERIFICATION evidence, not a restarted
trajectory. All physical workers are terminal. The fast direction has large
interior div(J) and radial normal-current jumps, with negligible wall normal
current. Volume-only entropy formulas do not reproduce the actual species rate.
No causal defect or repair has been established.

## Exact object to audit

The existing production species operator is K=(A-A^T)/2, where

    A = S^T W J_h · grad(N T),
    T = (N^T W N)^(-1) N^T W S,
    Cdot_s = -a_s K s_s / (e V),  a=(-beta, alpha).

S is the original fluid polynomial/Fourier reconstruction, N the compatible
continuous potential reconstruction, W physical cylindrical quadrature, and
V the actual fluid point volume. The default T is the physical L2 projection,
not the density-weighted historical point transfer. Its transpose has global
support and must not be interpreted as a compact face flux without proof.

At the stationary base v0=J0=0, the exact current-driven linear block is

    delta Cdot_s = -a_s K(delta J) s0_s / (e V).

No perturbed density/EOS derivative enters this block through J0. This exact
linear expression can be compared with the already archived finite-amplitude
centered species contribution. Retain both amplitudes/signs and conditioning
differences; do not replace any evolved state with this diagnostic expression.

## Independent full-vector weak decomposition

For every fluid test coefficient vector phi, let sh=S s, phih=S phi,
es=N T s-sh and ephi=N T phi-phih. Then

    phi^T K s = integral phih J·grad(sh)
                + 1/2 integral phih sh div(J)
                - 1/2 sum_cell_faces integral phih sh J·normal
                + 1/2 integral [phih J·grad(es) - sh J·grad(ephi)].

All integrals use actual one-sided reconstructed fields. Interior faces retain
both cell traces; periodic endpoints pair. Axis and wall contributions are
measured, not set to zero. The face sum uses cell outward normals. This identity
is a volume integration-by-parts decomposition, not an invented replacement
numerical flux. It must reconstruct the original K action before interpretation.

Implement its entire coefficient vector through independent tensor adjoints,
not only selected test functions. Use root-product fluid/field basis values and
derivatives, independently integrated L2 projection and its transpose, and the
audited explicit-cell normal traces. The fluid mass/volume is retained. Compare
production-quadrature projection versus independent q16/20 projection separately
so quadrature differences are not hidden as a changed physical operator.

The production scalar s_s is the point value C_s/n and sh=S(C_s/n). The previous
physical audit used (S C_s)/(S n). Record this representation difference and its
derivatives explicitly; do not silently identify the two definitions. For the
species-only linear block, physical entropy rate is S(delta Cdot_s)/(S n0).

Save separate advective, volume-charge, all-face, and transfer-residual vectors;
their native sum, physical reconstruction and independent component norms;
radial/axial face partitions; and the complete original/replayed input identity.
Use both the saved-error direction and independent smooth continuum mode.
Measure q16/20 convergence and require a fresh full replay without unit tests.

Only after this accounting should a candidate repair or field-refinement test
be selected. Any refined-field direction needs explicit B and current transfer
errors. Removing physical Hall/electron terms, projecting an accepted state,
adding filters/diffusion/floors or changing masses/bounds is not justified by
the current evidence. Full original formation/physical-validity gates remain.
