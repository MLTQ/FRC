# Kernel-restored balanced wall induction

Opt-in UniformWallInduction variant dispatching its inner solve to
BalancedWallCore. All full restored Ohm/Faraday, reconstruction and original
inner bounds, five outer attempts, physical masses and fixed preconditioner
are identical to the frozen implementation. The extra `attempt_trace` retains
all inner residual histories including rejected attempts. Numeric residual
components keep their original schema and acceptance meaning.

This wrapper is separate because the frozen StrongWallInduction explicitly
calls its original core and must remain reproducible. It changes only iterative
refinement; the wall Galerkin curl and separately measured raw-curl defect are
unchanged. Failure returns no accepted rate and does not mutate inputs.
