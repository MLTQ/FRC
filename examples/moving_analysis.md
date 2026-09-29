# Independent inspection of moving-plasma archives

`read_study` requires all declared cases complete, source/configuration matches
and verified lossless archive hashes. `compare` uses corresponding material
nodes on identical or nested doubled meshes and the coarser positive inertial
mass norm; it compares physical saved times, not frame indices alone. Relative
differences in small spurious equilibrium velocities are not physical accuracy.

`interior_flux_error` samples the actual Q1 interpolated psi at corresponding
material points. Its discrepancy measures the explicit difference between
nodal invariance and continuum frozen flux. It is never relabeled as physical
resistive reconnection. `expansion_pressure_error` integrates the exact expanding
gas pressure over actual final cells and compares cell averages; the exact
smooth formula is extended across small numerical surface-position errors.
`nonuniform_advection` independently prescribes a nonuniform two-dimensional
deformation and records refinement of the same interior-material-flux defect;
it is a component check, not an additional evolved plasma case.

Analysis reads but never repairs the accepted state. Its separate source hash
is retained by the report, avoiding retroactive changes to measured provenance.
