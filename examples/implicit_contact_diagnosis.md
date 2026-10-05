# Audit stopped inaccurate contact runs

Three Gauss candidates were stopped after substantial growth of deviations
from the exact translating contact. They did not complete the requested transit
and are not solver-validation successes. This audit hashes their last saved
accepted states and termination records, reconstructs the exact original
initial condition and independently measures both quadratures at the saved
checkpoint time. The saved checkpoint is not asserted to be the last accepted
state in memory when interrupted. No rejected trial is audited as accepted.

Reports full E, both temperatures, all field errors and physical energy using
the separately documented implicit_plasma_audit measure. Freezes source and
refuses overwrite. These observations do not distinguish spatial instability,
time aliasing, roundoff amplification or physical instability without further
independent controls. The source closure and failed pilot data remain preserved.
