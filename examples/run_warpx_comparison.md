# Bounded external-engine protocol

Runs 17 small verification cases sequentially on an explicitly selected GPU.
Requires the isolated optional WarpX environment and the two Python example
modules in `bench/examples`. No production solver or system package is changed.
Each case has a separate process, input record, log and output directory.

Nine base cases provide three meshes for both Hall helicities and Maxwell
vacuum. Eight controls halve and quarter each Hall time step, increase particle loading,
reduce amplitude, raise an inactive denominator floor, and halve the Maxwell
time step. The floor control is in uniformly positive plasma; it is not a
test of vacuum regularization. The Python report audits every output.

On restart the launcher skips directories containing a measurement record;
`backend_comparison_report` must still verify completion and source hashes.
Incomplete directories are deliberately not overwritten. Diagnose/preserve
failed trials, then choose a new run root or move that identified trial aside.
The shell launcher does not claim correctness from process exit alone.
