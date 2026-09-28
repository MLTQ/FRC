# Closed-flux equilibrium verification

Independent centered differences of continuous B and pressure check curl(B)×B
against grad(p), including field normalization and physical dimensions. The
power-series differential residual and boundary root are checked across the
declared aspect-ratio range, and Bz must reverse radially.

The default p0=1 Pa equilibrium permits affordable two-grid force,
finite-time flow and flux-drift checks, with conservation and discrete div B.
The p0=0.1 Pa coarse-grid case remains a reference failure with admissibility
limiting explicitly disabled; it must raise without any state repair. Separate
positivity tests verify that the enabled limiter admits this same input and
conserves cell budgets. Reference failures remain part of the evidence.
