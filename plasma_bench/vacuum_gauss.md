# Coupled Gauss evolution with free wall flux

VacuumCoordinates adds the physical outer psi row to the existing modal wall
coordinates. Bulk coordinates/scales, toroidal trace-null basis and real residual
norm are unchanged. There is no auxiliary current state: weak Ampere determines
all current coefficients from the magnetic field including the vacuum energy.

VacuumGauss inherits the existing nonautonomous DrivenModalGauss stage solve,
endpoint checks, source-time handling and accepted-stage source-work integration.
Its explicit field_type is VacuumMagnetic. No duplicated timestep algorithm,
operator splitting, step retry or changed nonlinear tolerance is introduced.

VacuumBlocks uses Fourier symbols extracted from the actual full geometry
matrices for its approximate uniform-state inverse. It includes constant physical
resistivity and current inertia as well as the existing ideal/Hall/fluid blocks.
The real nonlinear rates always use the full evolving state and source field.
All modes participate; no physical filtering or extra numerical damping occurs.

Five complete coupled runs and measured spatial/time/model comparisons are in
docs/VACUUM_CYLINDER_RESULTS.md. This bounded periodic coupling is complete;
finite coils/ends, moving free surfaces and formation remain open.
