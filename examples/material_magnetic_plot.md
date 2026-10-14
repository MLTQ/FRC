# Magnetized free-boundary scientific view

Reads complete r4/r8/r16 trajectories and reconstructs their actual final fields.
Produces one scientific PNG/PDF and quantitative refinement JSON; no new physical
runs, acceptance framework, tests or source edits. All field/current/velocity
comparisons use a common 32x64 material-label grid, not common Eulerian locations.
Weights are the finest sampled particle weights. The current is the weak-Ampere
projection including face terms; within-cell current is kept in run summaries.

Figure shows initial/final Bz and signed Jtheta, material flux contours, magnetic
energy/toroidal motion, and independent initial spherical-field/force errors.
Contours are the fixed material psi at vertices mapped to actual positions.
They are diagnostic flux surfaces of these initialized ideal data, not evidence
of coil-driven formation, resistive reconnection, stability or confinement.
No external field or coil exists in these source-free mathematical runs.

Initial and final panels share the same component color limits; no per-frame
color rescaling is used. Native cell values are volume-weighted for the display.
