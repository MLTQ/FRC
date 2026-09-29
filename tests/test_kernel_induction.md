# Kernel precision verification

An independently specified Er(r)/Ez(z) kernel must vanish exactly under the
difference-first physical curl; random fields agree with the sparse operator.
The complete original three-block magnetic/Ohm/Ampere equations are assembled
and directly solved separately from the magnetic elimination. Adding an electric
kernel to physical signals at amplitudes 1, 1e-4 and 1e-8 must retain those rates
within the declared input-roundoff amplification bound. There is no cutoff.

The small-signal solver keeps its strict relative residual; original full
equation backward errors are additionally checked. Full-E rate-relative errors
are descriptive conditioning diagnostics, not silently treated as tiny. Actual
contact evolution and source/work audits remain necessary beyond these tests.
