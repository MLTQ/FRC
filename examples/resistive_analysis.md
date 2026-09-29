# Independent comparisons of frozen resistive records

Requires completed, source/configuration-matched and hash-verified stationary
and moving archives. Comparisons retain the same physical diffusivity; changing
D is never described as a numerical accuracy improvement for the same problem.

Velocity uses corresponding material vertices and the coarser positive inertial
mass. Fine density/pressure are restricted by actual cell volume; integrated
cell heat is summed then divided by grouped volume. These compare the same
initial material sectors, whose finite geometric approximation also refines.

`material_samples` evaluates the actual Q1 field gradient and recovered current
at corresponding reference-sector Gauss points, including subdivided fine cells.
Common quadrature is composite over every fine subcell so fine-gradient jumps
are never integrated as one smooth coarse polynomial.
The two evolved physical coordinate maps can differ: this is a material-coordinate
state comparison, not interpolation onto an Eulerian Cartesian grid. Magnetic,
density, pressure and current-change norms subtract each mesh's initial state;
total-current error is also retained. Norms and absolute velocity/position
differences remain distinct, especially for small perturbations.

`restrict_cells` never repairs or modifies saved fields. Analysis source hashes
are separate from the frozen measurement source. An independent polynomial
reproduction check must validate field sampling and conservative restriction.
