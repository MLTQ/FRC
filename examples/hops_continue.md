# hops_continue.py

Continue the original user-authorized HOPS run without resetting its particles,
controller, source stream or budgets. No new physics or parameters are selected.

`main` freezes sources and starts a one-thread detached worker with nohup, its
own session and explicit log. A campaign file lock prevents duplicate workers.
The worker reconstructs legacy metadata using `iec_recover`, then replays the
last original 10 ns interval through the actual full 3D solver. Every saved
array must match bitwise (including dtype/signed zero), all physical diagnostics
must match exactly, and independent reconstructed RNG/controller states must
match the replay endpoint. Failure leaves the original interrupted.

After that gate, append from 2.52 to the original 5 µs target using a complete
checkpoint, preserving the original 253 frames and source archive. Each new
frame saves complete restart state atomically. Regular future restarts can use
that complete checkpoint directly. Geometry, physical code and config hashes
must still match. No invented emission, field repair or reset energy ledger.

DC/fixed controls follow sequentially with the original seed/source/parameters.
Completed cases are skipped; existing partial cases require full checkpoints.
Campaign/case STOP controls remain. Interrupted diagnostics are historical, not
a current error once continuation succeeds. Numerical/physical accuracy limits
remain those in HOPS.md; exact restart does not prove confinement or resonance.

CLI: `.venv/bin/python -m examples.hops_continue outputs/hops_resume_DATE
--run outputs/lab_runs/hops_v1_hops --campaign outputs/hops_v1`.
Read the recorded PID, `progress.json` and `worker.log`; detached launch alone
is not evidence that replay or physical continuation succeeded.
