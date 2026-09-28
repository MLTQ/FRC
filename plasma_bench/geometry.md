# Geometry
`Cylinder` describes chamber/plasma envelopes in meters. `Grid` creates uniform
nodal arrays including the axis and outer boundaries. Arrays are `(nr,nz)`.
Electromagnetics, solvers and visualization depend on this ordering. Dimensions
are finite positive values; at least five nodes permit second-order derivatives.
Region membership uses a roundoff-sized tolerance to preserve reflection symmetry
when decimal-valued endpoints lie on mesh nodes.
