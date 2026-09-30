# Actual physical derivatives and one-sided traces

Provides volume derivatives of the existing compatible electric/current and
regular-fluid fields, plus explicit left/right traces at field-cell faces.
Derivatives are of the declared polynomials/Fourier interpolants, not finite
differences of sampled pictures. Electric edge bases require second nodal
derivatives. Vector gradients use the physical cylindrical connection terms.
The degree-one edge derivative is exactly zero, without an out-of-range target.

`TensorSurface` supplies actual cylindrical surface weights. `current_maps`
evaluates a specified cell's polynomial at a face, so normal-current jumps are
retained. `fluid_map` similarly respects discontinuous radial fluid elements.
No epsilon displacement of a physical boundary or radius floor is used.
`RegularDifferential` differentiates the actual velocity quotient, including
the density dependence of axial velocity. Its `stress_rate` is the volume
gradient adjoint only; shared face tractions must be added by the caller.
