# Guarded complete fluid-refinement workflows

Fixed explicit dispatch names: `transfer` calls coupled_direction_transfer with
no additional input, `probe` calls coupled_fluid_probe with the verified transfer
directory and radial element count, and `audit` calls coupled_fluid_audit with
the completed probe directory. The latter modules are separate protocols and
must exist before their corresponding command is invoked.

Each selected module and this guard are explicit source seeds. Dynamic dispatch
does not silently omit the running scientific source, nor include future
unexecuted modules in an older identity. One-second own-process monitoring and
final extrema preserve 8 GiB RSS / 256 MiB free disk bounds. The first transfer
also requires a conservative 512 MiB full-study/replay payload reserve above
the disk minimum before starting. Current raw-array shape estimates fit within
that reserve; actual disk remains monitored throughout every workflow.

Optional `--original` performs the full saved-protocol replay verifier and raw
SHA comparisons. Reproduced physical RHS rejections remain rejections, not
accepted scientific cases. No time integration, retry, acceptance relaxation
or source replacement is introduced by the guard. Collect actual worker exits.
No unit tests; one heavy process at a time.
