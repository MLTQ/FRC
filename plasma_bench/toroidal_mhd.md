# Nonlinear compressible toroidal MHD core

`ToroidalMHD` integrates all three velocity/field components through six
conservative cell variables and corner psi. Reuses the validated mesh,
parameters, poloidal resistivity and scalar trace limiters; the existing
four-variable MHD solver is unchanged. Full toroidal stress, energy and fast
speeds enter fluxes, source integration and UCT. Radial wall is impermeable,
slip and conducting; z periodic; every cell contains fluid. No vacuum or
formation interface is inferred. Hall and finite electron mass are not yet
nonlinearly coupled.

SSPRK2 validates both trial stages before atomically committing state, time,
heat and limiter accounting. No pressure floor or energy reset exists. CFL
combines full fast waves and a conservative bound on explicit resistive
poloidal/toroidal diffusion. `advance` supports a declared independent time cap
and preserves accepted partial state on a step-budget exception. Diagnostics
report angular momentum and toroidal flux separately from total energy and
positive physical heat. Total energy includes Joule conversion once, not twice.

The first/second-order reconstruction and Rusanov jumps are explicit numerical
stabilization. Positive face traces do not establish full-update robustness.
The compatible radial current mass recovers regular polynomial diffusion at
the axis without sacrificing positive work. Remaining smooth/shock errors must
be measured. Equations and verification: TOROIDAL_MHD_NUMERICS and
TOROIDAL_MHD_VERIFICATION.
