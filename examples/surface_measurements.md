# Resolving the material surface layer

Records twenty moving runs at unchanged D=0.1 m²/s, gamma=5/3 and rho=1 kg/m³:
source meshes nr=8/16/32 with grading 1/1.5/2, unforced resistive FRC relaxation,
and independent angular, outer-time/CFL, diffusion-tolerance, magnetic-matching
and quadrature controls. Physical inputs match the preceding resistive fixture.
Grading changes both material and vacuum shell spacing near the interface.

Each case uses the shared complete resistive recorder with this study's full
fingerprint, including actual case and reference helpers. Reuse checks fields
as well as complete configuration and sources. Disjoint names can run in
separate processes; a partial file is progress evidence, not a live-job handle.
`--references` separately records thirteen independent gas/decay controls.
The historical uniform-grid study remains archived at its original source.

Neither grading nor positive heating resolves missing surface thermodynamics
automatically. Report all field, current, pressure/density and velocity profiles,
including changing flow sign at the surface, and quantify remaining uncertainty.
