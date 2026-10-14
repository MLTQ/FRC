# Coupled refinement boundary at the display request

**Superseded 14 October 2026:** F32 physical audit and the complete comparison,
including their full replays, are now finished. Read
[COUPLED_FLUID_REFINEMENT_RESULTS.md](COUPLED_FLUID_REFINEMENT_RESULTS.md).
The text below preserves the earlier display-request boundary.

The user requested an updated 3D field/current display while this study was in
progress. This note preserves its actual boundary. No additional physical
calculation was launched for the display. Read COUPLED_FLUID_REFINEMENT_DESIGN.md
and the new modules' companion documentation before continuation.

## Completed evidence

- `coupled_direction_transfer` and its complete replay subdivide both archived
  conservative directions onto F16/64 and F32/64 with original B32/64 and
  normalization. Raw contexts, transfer comparisons and replay are preserved.
- `coupled_fluid_probe_F16` and `_F32`, their complete fresh replays and verifiers
  contain all signed full-RHS calls at three amplitudes for two directions.
  Each probe has 26 raw files; every scientific replay comparison passes.
- `coupled_fluid_audit_F16` and its complete replay/verifier contain independent
  physical reconstruction, F8 comparisons and instantaneous physical integral
  rates. Both raw files and the complete metrics match byte for byte.

The display assembler checks these raw files, summary hashes and source
identities. New used scientific files are frozen: `subdivide_fluid_direction`,
`coupled_direction_transfer`, `guarded_coupled_refinement`, `coupled_fluid_probe`,
`coupled_rate_physical`, and `coupled_fluid_audit`. Earlier physical sources,
states, residual bounds and the user's poster remain unchanged.

Prior working context recorded actual zero exits for the transfer pair, both
probe pairs and F16 audit primary. The F16 audit replay session 30234 was no
longer available when inspected for this display request. Its saved guard says
terminal, resource-compliant and complete-protocol replay verified; its verifier
and raw hashes independently pass. Do not claim a newly observed process exit
for that expired session, or rerun it merely to recover a process handle.

## Next scientific action

F32 physical audit has **not** been launched. The prepared
`coupled_refinement_comparison.py` has **not** been executed. Finish the F32
audit and full replay, then run the complete physical F8/F16/F32 comparison
under a new dedicated guarded wrapper with a full replay. Do not alter the
already-used guard to add a dispatch case. One heavy worker at a time; retain
8 GiB RSS / 256 MiB free-disk bounds and original equation tolerances.

```sh
PYTHONDONTWRITEBYTECODE=1 VECLIB_MAXIMUM_THREADS=1 \
MPLCONFIGDIR=/private/tmp/plasma-bench-mpl MPLBACKEND=Agg \
.venv/bin/python -m examples.guarded_coupled_refinement audit \
  outputs/coupled_fluid_audit_F32 \
  outputs/coupled_refinement_diagnostics/F32_audit_control \
  outputs/coupled_fluid_probe_F32
```

Then inspect before launching the replay into a distinct destination. The
comparison must include original F8 physical energy rates and actual field
differences, not differences between error norms. Scalar-source refinement and
changed coordinate projections do not demonstrate a coupled cure or physical
instability. No newer accepted trajectory exists; step 17 at 28.3333 µs is still
the correct field display. Full signed/time/arithmetic trajectories, local
budgets, vacuum, transport, physical coils/ends and formation gates remain open.

The local repository has no remote. Export beads locally; do not publish data.
The display request does not close or resume the overarching goal by itself.
At this boundary the partial scientific continuation and new display are saved
in the working tree and remain uncommitted; the previous commit is 887faca.
Beads have been exported locally. Preserve these files and the separately
modified user poster when preparing the next scientific checkpoint.
