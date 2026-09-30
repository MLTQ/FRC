# Experimental radial pressure pairing

Separate diagnostic subclass preserves the frozen finite-volume baseline.
Replaces its radial pressure/pressure-hoop and centrifugal source by the physical
thermal/kinetic-work adjoint. Conservative convective transport, conservative
axial/angular stress, continuity, both entropies, energy account and RK4 remain
unchanged. The temporary duplicate derivative orchestration keeps all baseline
evidence replayable; it should be consolidated only after an adoption decision.

`radial_kinetic_work` differentiates the exact convective finite-volume flux
power with respect to the radial particle flux, holding transported specific
momentum traces fixed. It includes the actual kinetic density derivative.
Subtracting this kinetic work with the same metric supplies the cylindrical
centrifugal force and cancels radial advection energy defects. The separate
pressure-only trial failed rotating balance and is preserved in outputs.

No exact full energy theorem follows: axial pressure and convective work defects
remain measured, and the radial adjoint still needs independent force/rotating
equilibrium and long-trajectory validation. This is a derived semidiscrete work
pair, not a post-step energy repair. No diffusion or floor is added.
Measured original-equation rates expose a nonconverging relative radial force
error in the first annulus, despite small rotating-equilibrium residuals and
improved global norms. This variant remains experimental; see RADIAL_WORK_NUMERICS.
