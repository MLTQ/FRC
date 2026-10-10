# Construct only the metric used by the split-inverse solver

WallInductionCore constructs a global magnetic-Hessian LU. Its subclass
UniformWallInduction immediately replaces that factor with UniformWallMetric.
DirectMetricWallInduction copies the exact original validation, active indices,
null-space basis, normalization, lift and dual-lift expressions, then constructs
the unchanged UniformWallMetric directly. Only the discarded LU is omitted.

Every original split-inverse solve, physical mass, actual density, operator,
kernel restoration, iteration/refinement budget, residual and acceptance bound
is inherited unchanged. No accepted numerical coordinate or field is projected
or corrected. The uniform metric remains an approximate inverse only; actual
nonuniform density remains in the complete operator. No diffusion is added.

DirectMetricWallPlasma retains identical original fluid/current/species terms
and cache-copy rebinding to the actual current model. direct_metric constructs
a separate assembly; frozen original systems and sources stay unchanged.

Verification must use full original-versus-new signed workflows, requiring
bitwise complete inputs, rates, solver diagnostics and independently audited
fields, followed by fresh full new-workflow replay. Omitting the unused work
is not assumed to improve runtime until measured. No unit tests are required
under the user's current preference. This bounded model does not add vacuum.
