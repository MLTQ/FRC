# Complete pulse with physical wall-current evolution

Runs the same24µs pulse and relaxation as driven_flux_pulse: R35mm,L240mm,
Ar+,n1e21/m³,Te5/Ti1eV,peak mean axial field0.1mT for12µs. Only the physical
current space and associated dynamics change. The new state carries the
azimuthal wall-current coefficient, evolved with physical finite inertia.
Every original fluid, Hall, pressure, transport, exchange and Joule term remains.

Use fixed eta or actual evolving Spitzer resistance. The fixed-eta continuum
Bessel reference and its nonzero wall current remain independent of production.
Spitzer differences from it are control/model differences, not numerical errors.
Independent polynomial reconstruction includes the additional wall nodal value,
with no current field substituted into dynamics or heating. Native and physical
energy include its current kinetic energy; source work uses accepted Gauss stages.

Every endpoint saves y alongside fluid and B, plus the entire stage changes.
Full Ohm, original induction, stage and endpoint bounds remain. Output also
records actual and reference wall current, independently reconstructed wall B,
wall Poynting power, and its difference from discrete source work. These are
physical measurements of this solver change, not additional acceptance layers.

Own4GiB RSS/256MiB disk checks, source freezing and terminal records use the
existing driver pattern. No unit tests, new framework, numerical diffusion,
floors, heat repair or atmosphere. Filled periodic rigid cylinder only; no
free plasma edge, solved vacuum, finite coils, real ends or FRC formation.

Measured full24µs160/320-step runs and their remaining spatial/current limits
are recorded in docs/WALL_CURRENT_RESULTS.md. The finer spatial run is still
active; this is not yet an accepted formation or vacuum-boundary model.
