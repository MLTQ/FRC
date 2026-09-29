# Explicit unsplit finite-mass time integration

Classical RK4 advances conservative mass/momenta/total energy, psi/g and electron
U through the same four complete coupled RHS evaluations. It never calls a
fixed-n substep or reinitializes total energy. Every stage and final species
partition must validate before state, time or accepted budgets commit. Failed
trials raise without a smaller hidden retry, pressure floor or state repair.

The caller must prescribe max_timestep and establish resolution/stability with
independent time refinement. This cap is not an automatic dispersive CFL proof;
explicit electron inertia/Hall and drift dynamics may require very small steps.
RK4 has fourth-order smooth consistency, but full measured convergence and
stability are required. A future implicit method may improve efficiency after
this coupled reference is verified. No physical diffusion is added here.

Diagnostics expose conservative integrals, external/boundary energy input using
the actual RK quadrature, positive species states and linear residuals. Global
energy conservation does not prove local temperature or entropy accuracy. Source
callbacks and initial_time are explicit; unforced runs omit the callback.
Nodal electron and cell ion log-entropy diagnostics use their stated discrete
representations and a fixed density/energy gauge; they are not an entropy
inequality or an independent physical quadrature. Partial accepted state remains
available if an advance step budget is exhausted. Stage2 is not enabled and
the physical plasma-vacuum surface is not represented by this filled cylinder.
