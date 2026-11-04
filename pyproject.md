# pyproject.toml

Package metadata and ordinary CPU dependencies for plasma-bench. The
`precision` extra selects existing arbitrary-precision work; `dev` retains
historical pytest tooling, which is not run under current project instructions.

The optional `metal` extra pins MLX 0.32.4 on Apple-silicon macOS. Its pinned
mlx-metal dependency supplies the actual GPU runtime. GPU availability is
checked only when explicitly selecting the experimental backend; normal
CPU imports/runs do not require MLX. Installed wheel requirements still
apply (the measured Python 3.13 Metal wheel requires macOS 26).

The optional `electrostatic` extra adds Gmsh 4.15.2 for conforming 3D cage meshes
and Numba for tetrahedral particle traversal. It does not change the FRC solver.
Generated mesh arrays and ordinary SciPy float64 sparse Poisson solves remain
local; there is no cloud runtime or GPU claim for this initial electrostatic path.
