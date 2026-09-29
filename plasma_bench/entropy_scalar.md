# Physical scalar metric for specific entropy

`EntropyScalar` samples continuous nodal degree-3/5 scalars with the same Nx,Nz
as the compatible fields. It takes the **existing** immutable conservative
density map: n=R n_cell. The physical metric is M_n=N.T W n N, with cylindrical
quadrature W. No lumped dual volumes or alternate nodal density are substituted.

Tensor actions avoid a large multidimensional sampling matrix. The unweighted
physical scalar Gram matrix preconditions a CG solve of the full variable-n
metric; normalizing by volume-mean n changes only numerical scaling. Every solve
checks its original load residual and returns iteration count/tolerance evidence.
Failure is atomic. Constant projection anchors are an exact linear decomposition,
not an equilibrium fit. A zero load returns an exactly zero rate.

`gradient` uses the compatible incidence before field sampling. Constant scalars
therefore have exactly zero gradient, with no axis division of a scalar. `adjoint`
is N.T W including volume. `mass_action` retains physical n units. `project` is
the number-weighted L2 projection, not point interpolation or a cell average.

This scalar space alone supplies no continuity law, thermal model, transport
limiter, shock entropy inequality, vacuum carriers or time integration. Positive
quadrature density is inherited from `CompatibleDensity`; no second density is
evolved. The local overlapping polynomial basis has no generic stability theorem.
