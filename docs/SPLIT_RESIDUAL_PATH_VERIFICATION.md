# Residual arithmetic and component response continuation

The full FRC goal remains open. No unit tests were created or run. These are
filled-periodic verification fixtures, not selected laboratory conditions.
See SPLIT_RESIDUAL_PATH_DESIGN.md for equations and unchanged assumptions.

## Exact input observation

The F64/256, B256/256 rejected minus input is reconstructed bitwise. The
observation-only current split-inverse core reproduces the original first-
attempt failure exactly, including final operator 3.2861098689e-11, dual
Faraday 5.8882016542e-11 and physical Faraday 1.5875361928e-11 under the 1e-11
bound. Full Ohm/Ampere backward errors remain 5.3805e-16 / 2.7490e-16.
No failed plasma rate is accepted and no solver arithmetic or bound changed.

The measured opposing residual-vector gap is 6.46465e-11 at the initial
iterate and 6.87234e-11 after the third correction, normalized by the original
right-hand side. The electric curl-distribution term alone leaves only
0.4245% of the final vector gap unexplained. H/T metric association is
6.22e-15 and inertia-path discrepancy 2.89e-13 in the same normalization;
neither dominates this rejection. All seven measured terms together leave
5.41e-17 of the right-hand-side norm, or 7.87e-7 of the gap.

An independent exact Decimal audit checked all 196608 electric coefficients
at each of four iterates (786432 coefficients total), with zero mismatches
against the recorded two-sum low terms. Precision is 2200 decimal digits;
all stored binary64 values are converted exactly. The complete audit and its
fresh replay agree bitwise. Nevertheless, propagating the lost electric-addition
bits alone leaves 66.21% of the final residual-gap vector unexplained. Matrix
curl distributivity must remain separate from addition loss. A term norm alone
is not evidence that it predicts the observed vector discrepancy.

The full observation replay was interrupted by its own disk guard at
197836800 bytes free during CG, after the original observation had completed.
Its source, partial summary and traceback are preserved. This is a resource
interruption, not a numerical rejection or a completed observation replay.
The fresh standalone replay subsequently completed and matches every saved
array, iterate observation and non-runtime metric bitwise. No unrelated files were
deleted. The original observation occupies about 88 MiB; exact arithmetic
audits about 2.3 MiB each. The much larger host disk-headroom change is not
explained by those output sizes alone; no system-wide cause is asserted.

## Physical component reconstruction

Both full physical reconstruction workflows and their fresh replay are
complete. All original physical-rate, E, J_t, reference and quadrature-weight
fingerprints match. Whole-volume relative L2 errors and actual differences
are normalized by each component's own independent continuum reference:

| Response | Error at epsilon 1e-4 | Error at epsilon 1e-5 | Physical difference |
|---|---:|---:|---:|
| Radial current rate | 0.05969389 | 0.06672077 | 0.03014823 |
| Azimuthal current rate | 0.0001904658 | 0.0001904658 | 4.52804e-10 |
| Axial current rate | 0.01187001 | 0.01193623 | 0.001262864 |
| Full current-rate vector | 0.0008789930 | 0.0008854397 | 0.0001074501 |
| Radial velocity rate | 0.001831426 | 0.001833905 | 0.0001078761 |
| Toroidal magnetic rate | 0.005475095 | 0.005502829 | 0.0005524026 |
| Full electric field | 4.568074e-6 | 4.568896e-6 | 7.148893e-8 |

The total current norm therefore hides a 5.97–6.67% radial-component error and
about 1.19% axial-component error. The earlier sub-0.1% vector error remains
correct in its stated norm but does not certify individual components.
All other fluid, species-entropy, magnetic and electric components, physical
maxima and per-row squared integrals are retained in the full JSON/NPZ records.

At quadrature 20, the fixed inner/outer 5%-radius bands and remaining interior
contain the following fractions of squared amplitude difference:

| Response | Axis | Interior | Wall |
|---|---:|---:|---:|
| Radial current rate | 0.210% | 88.674% | 11.115% |
| Axial current rate | 0.300% | 88.188% | 11.513% |
| Radial velocity rate | 0.492% | 90.277% | 9.231% |
| Toroidal magnetic rate | 0.261% | 89.711% | 10.028% |

The radial-current amplitude sensitivity is not confined to the wall. Radial
current error is about 5.9–6.0% in each band at epsilon 1e-4. In contrast,
58.08% of squared axial-current error and 72.05% of squared radial-velocity
error lie in the outer 5% radius at that amplitude. Those are error fractions,
not amplitude-difference fractions, and the two must remain distinct.

Whole-domain q16/20 relative norm changes are at most 3.44e-15 absolute. Bands
mask the original Gauss nodes; cells are not subdivided at band boundaries.
The largest band-fraction change is 0.001274, for the precision-scale density
amplitude difference (whole-domain norm 9.56e-14). The radial-current wall
sensitivity fraction changes by 0.000190. These quadrature differences remain
visible rather than implying exact band integrals or continuous maximum bounds.

## Evidence and next work

- Original exact observation: outputs/split_residual_probe.
- Independent exact arithmetic and replay: outputs/split_residual_arithmetic_audit,
  its _replay and _verification directories.
- Complete component audit and replay: outputs/component_amplitude_localization,
  its _replay and _verification directories.
- Scientific figure: outputs/split_precision_report/verification.png, visually
  inspected; vector PDF alongside it. The PDF was not separately inspected.
- Interrupted original replay and resource record: outputs/split_residual_diagnostics.

All workers are terminal and collected. The fresh standalone full observation
replay and exact verification completed; its disk guard did not trigger. Earlier
interrupted data remains preserved. Next compare the projected curl applied
as an assembled matrix versus curl followed by projection, retaining all
original residual checks and reporting arithmetic discrepancies explicitly.
Do not infer that either factorization or addition compensation fixes the
complete plasma solve. A causal connection between this residual roundoff and
the physical current-component errors has not been established. Also separate component error by axial/radial modes
or resolution before adopting the previous total-current accuracy claim.

Full signed trajectory convergence, local budgets, moving exact vacuum,
transport/ionization, real coil/end conditions, converged prescribed-waveform
formation and the experiment UI remain open. No numerical diffusion, artificial
mass, atmosphere, field filter, heat/energy repair, remote or publication was
introduced. No unit tests were added or run. All earlier Python sources and
archived calculation evidence remain unchanged.

## Final preservation checkpoint

Three complete verification groups cover the full captured-state observation,
independent exact arithmetic audit and physical component reconstruction.
All saved arrays and non-runtime metrics replay bitwise. Twelve new source
archives and all 288 preceding archives are verified; 22 older source closures
still match. All six new Python modules have companions; no older Python source
changed. Final 134-file source union:
150dc2f3e7614eed5da58d314b6ff3e8cf804b641c1f9fec9f91b738e4712745.
The interrupted replay, guard record, original scientific evidence and user's
poster are preserved. Full inventory: outputs/split_residual_diagnostics/checkpoint.json.
All 12 acceptance beads and the overarching FRC goal remain open.
