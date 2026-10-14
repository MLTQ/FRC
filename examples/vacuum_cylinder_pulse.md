# Full plasma/vacuum pulse and optional axial mode

Runs the full 24 microsecond pulse/relaxation with R=35 mm, L=240 mm, ideal
cylindrical source radius C=75 mm, Ar+, n=1e21/m³, Te/Ti=5/1 eV and gamma=5/3.
The source field peaks at 0.1 mT for a 12 microsecond sin² pulse. Plasma flux is
free, exterior matter is absent, and the source sheet has prescribed current.
Select fixed or evolving Spitzer resistance. Optional seed_amplitude adds the
independent lowest nonconstant axial vacuum-matched decay mode to the initial
field, exercising nonzero exterior induced field and energy during the pulse.

Uses the actual coupled VacuumPlasma/VacuumGauss solver. Saves every accepted
state, full stages, source-work history and frozen source closure. Independent
root-polynomial reconstruction measures plasma field/current/heating and mass.
Exterior field is compared against the separate continuum K-Bessel solution;
physical energy sums retain vacuum induced and source-background energy. Full
native energy also exposes the source/induced cross term. Field/current errors
are against a linear fixed-resistance baseline, so nonlinear/Spitzer differences
must not automatically be called numerical errors.

Original own 4 GiB RSS / 256 MiB free-disk guards and solver acceptance bounds
remain. No unit tests, replay framework, artificial atmosphere, numerical
diffusion or heat repair. Fixed footprint and periodic axial boundaries remain;
this run does not establish real-end behavior, formation or confinement.

Radial and matching fluid/field axial counts are explicit nr/nz options (8,16,32).
The initial uniform 8/16 run predates the nz option; its frozen source is retained.
The option changes spatial resolution only, not equations or acceptance bounds.

Completed runs, source closures, error norms and measured limitations are in
docs/VACUUM_CYLINDER_RESULTS.md. The axial seed is initialized as a mathematical
mode; any initial closed flux is not evidence of formation.
