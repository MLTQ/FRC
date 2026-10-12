# Serial isolated full physical controls

`run(destination, plan_path)` executes explicit guarded_physical_controls jobs
one at a time, each through the same Python interpreter in a fresh child
process. Plan entries give name, physical/control destinations, nr, nz, epsilon,
tolerance and optional original directory for complete replay. The source and
exact plan hash freeze before execution. Existing campaign paths reject.

Each child retains the original 8 GiB/256 MiB guards and exact physical protocol.
The parent saves actual child PID, status, return code, log hash and terminal
control record. A child resource/error exit stops subsequent jobs. A terminal
strict numerical rejection is explicitly recorded and can be fully replayed;
it is not an accepted physical rate. Parent interruption signals only its own
live child, collects it, and preserves all partial evidence. No shell is used.

`share_verified_replay` runs only after a complete independent replay, validates
verification/summary/raw hashes, and atomically shares byte-identical NPZs.
Every original path and byte remains, with a full manifest. Shared paths are
immutable; never rerun into an existing output directory or write either path
in place. No unreplayed, differing or interrupted data are removed.

The campaign is execution bookkeeping, not an additional physical validation
or a substitute for its individual simulation replays and independent audits.
Do not run another solver/audit concurrently. No unit tests, changed equations,
relaxed numerical bound, new atmosphere/filter/floor or laboratory prescription.
All sources and completed evidence stay frozen for subsequent continuation.
