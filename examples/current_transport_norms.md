# Independent transport norms and physical power

Evaluates all three saved electric transport fields against analytic original
covariant terms. Drift force is reconstructed from the actual velocity-state
rate (including density quotients), not merely compared at the solver nodes.
Reports global, first/last fluid-element and maximum errors, with separate
radial/azimuthal/axial force components and omitted-face negative controls.

Independent physical integrals measure density exchange, current/magnetic
energy rates, current-electric work, true strong continuity and div(J). The
full reconstructed fluid energy derivative is reported separately from the
grid energy identity; finite-resolution differences must remain visible.
Internal energy is n*e*exp(C/n)*(n/n0)**(gamma-1)/(gamma-1), using the
original entropy convention and the declared reference density.
No power is subtracted to force conservation. Higher quadrature measures
integration error, not physical discretization error.
