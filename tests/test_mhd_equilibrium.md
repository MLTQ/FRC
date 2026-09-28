# Closed-flux equilibrium verification

Independent centered differences of continuous B and pressure check curl(B)×B
against grad(p), including field normalization and physical dimensions. The
power-series differential residual and boundary root are checked across the
declared aspect-ratio range, and Bz must reverse radially.

The default p0=1 Pa equilibrium permits affordable two-grid force,
finite-time flow and flux-drift checks, with conservation and discrete div B.
The p0=0.1 Pa coarse-grid case is a known reconstruction failure and must
raise without any state repair. The full report retains finer-grid drift and
failed attempts for that more demanding case; it is not omitted from evidence.
