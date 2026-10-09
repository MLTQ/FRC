# Independent high-radial physical-field check

Compares the selected radial16,N128 reference current with the original cylindrical
curl(B)/mu0 computed by an independent fourth-order centered finite difference
at three physical locations. Checks all three components, including cylindrical
Btheta/r. It does not compare the loader to its own internal coefficients or
fit a production eigenvector. Axial eigenfunction convergence is a separate
archived numerical study.

The original N64 selection failed this unchanged tolerance with0.003251 error
in Jr; independently refined128/256 references pass. This is a reference
resolution requirement, not a loosened test or modification of the old archive.
