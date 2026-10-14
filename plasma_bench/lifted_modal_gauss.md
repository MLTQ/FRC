# Add magnetic states in admissible coordinates

`LiftedModalGauss` retains the original complete GuardedModalGauss equations,
physical scaling, Gauss coefficients, inverse, nonlinear proposals, line search,
stage/endpoint residual checks and strict induction/wall bounds. The frozen
original implementation is preserved; this opt-in variant changes how magnetic
stage and endpoint states are constructed in floating-point arithmetic.

The original method forms g_stage=g_initial+Q delta_y, where Q spans the
toroidal wall-trace nullspace. Nearly cancelling admissible vectors can lose
relative wall accuracy: the failed Spitzer8/160 endpoint had normalized trace
5.87858e-14 against the unchanged5.68434e-14 limit, with coefficient norm
4.14217e-16 and absolute trace norm1.75477e-28. Its complete probe is retained
in outputs/spitzer_bessel_decay_r8_t160_boundary_probe. This identifies the
immediate cancellation rejection, not the separate long-contact accuracy bug.

The variant stores independent magnetic coordinates y_initial and constructs
g_stage=Q(y_initial+delta_y) with a single lift. Poloidal coefficients follow
the same coordinate addition. No already constructed field is projected,
clipped or overwritten to pass a check. All physical rate and accepted-state
differences are still measured against the actual input field by the original
residual formulas. The coordinate norm, basis and acceptance bounds are fixed.

On a fresh input, admissible magnetic coordinates are encoded once. Thereafter
the exact accepted coordinate array is cached alongside its physical field;
it is reused only when the supplied field matches byte for byte. Failure does
not commit the candidate coordinate array. Arbitrary different inputs require
a fresh encoding and remain subject to the original checks. Fluid arithmetic
is unchanged. This is a representation choice, not a new physical model,
numerical damping term or energy correction.

The step body is intentionally a frozen-source variant of GuardedModalGauss
because old archived trajectories must keep their implementation. Only the
coordinate cache and magnetic construction differ. Use the full Spitzer decay
example through examples/lifted_spitzer_decay.py for end-to-end evaluation;
no unit tests or new verification framework are introduced.

Both full12us Spitzer80/160-step runs complete at the original bounds. Energy
drift decreases15.9858-fold; old/new80-step final vector-B difference is
1.90935e-10 of initial Bz in independent physical L2 norm. This supports this
arithmetic repair, not a cure for the separate long-contact growth. Full
results and preserved original rejection are in SPITZER_DECAY_RESULTS.md.
