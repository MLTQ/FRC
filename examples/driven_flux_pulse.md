# Driven collisional cylinder and relaxation

Evolves one12us sin² enclosed-flux pulse followed by12us relaxation in the
complete moving finite-mass model. R0.035m,L0.24m,n1e21/m³,Ar+,Te5/Ti1eV,
both gamma5/3, initially zero velocity/current/field. Peak mean axial field
is0.1mT. This is a mathematical drive case, not selected laboratory geometry
or preionization. The prescribed wall loop voltage follows Faraday exactly.

`fixed` uses eta3.8318200089006576e-5ohm m and Coulomb heat exchange;
`spitzer` uses actual evolving weak-field Z1 resistance with the same heat
exchange. Electron Joule heat, motion and every original two-fluid term
remain. The independent constant-eta reference includes current inertia and
the correct enclosed-flux boundary. In spitzer mode its differences are a
comparison with a different physical model, not numerical error estimates.

Read flux_diffusion_reference.md for the continuum formula and the analytic
polynomial that retains nonzero wall current. Selected-phase512/1024-mode
differences measure truncation. Production current has the old zero wall
trace, an explicit boundary approximation whose physical error is measured.
No model or accepted state is changed to match the reference.

DrivenModalGauss uses actual source times and fields at both Gauss nodes and
the endpoint. Accepted stage powers supply integrated source work. Independent
order20 reconstruction measures full B/J, fluid state, energy and particle
number; both physical and native energy minus source work are scaled by the
peak harmonic-field energy V A²/(2mu0), not the large thermal background.
All accepted states/stages, source closure and original failure evidence are
saved. Source-time arithmetic is carried forward from each accepted endpoint.

The exact mean axial flux is an imposed constraint, not evidence of numerical
accuracy. Differences from its physical integral are recorded separately.
The scalar resistance closure remains conditional; current/thermal drift and
magnetization diagnostics are retained in spitzer mode. Thermal conduction,
thermal forces, viscosity, neutrals and ionization are absent. The filled
periodic impermeable cylinder has no solved vacuum, finite coils, real ends
or free plasma edge. Closed contours or field reversal here are not formation.

nr8/16/32 is supported with half as many fluid radial elements, axial16 and
degree5. Original residual/boundary/positivity bounds remain; no numerical
damping, atmosphere, floor, retry or energy repair. Own4GiB RSS and256MiB disk
checks are discrete. No unit tests or additional verification framework.

Three full driven runs are recorded in docs/DRIVEN_FLUX_RESULTS.md. Fixed-eta
Bz improves on8→16 cells, but space-time current error remains28.06%→21.14%
and heating7.87%→4.47% low. The evolving-eta case is a physical control
comparison. The wall-current approximation remains an unresolved solver issue.
