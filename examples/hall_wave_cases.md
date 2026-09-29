# Hall cylinder measurement contracts

Defines 49 complete configurations: 24 mesh cases (three explicit models,
both helicities, nr=8/16/32/64), six collisionless controls, nine time controls,
two axial sample controls, two quadrature controls and six higher radial/axial
mode cases. MHD/Hall cases have R=1 m; finite-mass cases have R=.01 m to resolve
electron-inertia dispersion. These are independently named mathematical fixtures,
not alternative lab parameters or replacements for the original argon target.

Each case initializes by nodal interpolation of analytic potentials and fields,
then enforces their exact zero traces. It does not select a numerical eigenmode
or adjust the analytic frequency to improve agreement. Final time corresponds
to two radians of the continuum positive-frequency mode; all SI coefficients,
mode numbers, quadrature and time controls are explicitly recorded.

Three additional slow-branch time controls (2048/4096/8192 steps) resolve
secondary fast modes exposed by the initial 128–1024-step measurements.
