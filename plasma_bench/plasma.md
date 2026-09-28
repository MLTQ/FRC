# Plasma solver interface
`PlasmaSolver.step(external_psi, dt)` decouples prescribed electromagnetics from
plasma dynamics. `StationaryDiffusion` solves the mass/stiffness differential-
algebraic system documented in PHYSICS.md, with exactly zero mass in vacuum.
`DiffusionStep` returns induced flux, interval electric field/current, induced
energy and step heat/work/numerical loss/residual. Its energy identity is exact
for the discrete problem; it does not imply temporal convergence. Cached LU is
reused at fixed dt. Initial vacuum constraints and axis gauge are checked.
`ResistiveMHD` gates the original formation experiment until verified dynamic
plasma/vacuum coupling exists. The separate `mhd.AxisymmetricMHD` kernel has
conservative dynamics in an explicitly bounded verification cylinder; its
boundaries cannot be silently substituted for the quartz/vacuum experiment.
