# Fixed radial16 rates with split mass inverse

Separate opt-in copy of frozen balanced_radial_rates, changing only the full
plasma assembly to SplitInverseWallPlasma. Same physical N128 independent
reference, signed inputs, all original equation/residual checks, quadrature
16/20, full field/current/face and initial representation diagnostics. Full
coefficient/metric replay is required. Explicit failures remain terminal
failures. No old accepted solver, archive or complete trajectory adopts this
proposal. The electric reconstruction changes floating-point evaluation of a
linear inverse, not the physics, tolerances, spatial approximation or diffusion.
