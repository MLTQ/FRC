# Physical L2 transfer for species potentials

`SpeciesProjection` maps the regular scalar polynomial/Fourier field S f to a
compatible field potential N phi by physical unweighted L2 projection:
T=(N^T W N)^-1 N^T W S. Tensor-product Gram matrices make this a pair of small
direct 1-D solves; `adjoint` is the exact transpose of those saved maps.
Inverse Gram matrices give global radial support; axial Fourier interpolation
also has global support. Graph conservation is not a geometric face-flux proof.

This common integration is essential for consistency of the skew entropy
operator. A point-interpolation map's transpose does not distribute the weak
divergence consistently back onto the fluid metric. The latter is preserved
as `SpeciesPotential` for the failed-proposal control.

T preserves constants to solve roundoff and uses both element traces through
volume integration. It is a continuous common potential, not a new density,
temperature, entropy state or positivity-preserving interpolation. The physical
specific entropy is represented here by S(C/n); the independent continuous
diagnostic still reconstructs (S C)/(S n) and measures the difference. No extra
mass, filter, diffusion or global residual correction is introduced.

Dependencies: RegularSpecies uses `apply` and `adjoint` with actual regular
scalar shapes and full compatible nodal potential shapes including axis/wall.
