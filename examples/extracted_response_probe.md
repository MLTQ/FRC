# Actual spatial rate of the extracted fast response

Uses hash-checked signed trajectory snapshots at 40%,60%,80%,100% of duration.
Subtracts exp(lambda*40% duration) times the earlier response from each later
response to cancel the independent slow reference. This is diagnostic linear
combination of saved conservative coefficients; no evolution state is changed.
The two windows are separated by 20% duration. Physical component projections
give inferred growth and the remaining shape difference between windows.

Normalizes the late difference in fixed original physical coordinates and
evaluates the actual full strict stationary-contact RHS with both signs at
epsilon1e-4 and1e-5. Reconstructs its tangent rate independently at quadrature
16/20, comparing projected growth, residual not parallel to the direction and
amplitude dependence. This tests the spatial RHS without a time integrator;
it is not a full eigenanalysis or independent continuum verification.

Saves original coefficient directions, both rates, all localization arrays,
physical fingerprints, input hashes, strict residuals and source/runtime.
Requires a complete paired physical audit; full replay is separate. A nearly
eigenlike extracted direction does not identify its physical/numerical cause.
