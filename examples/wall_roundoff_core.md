# Instrumented unchanged balanced iteration

Exact opt-in diagnostic copy of BalancedWallCore solve arithmetic and bounds,
adding observations only. Retains each proposed coordinate vector, operator
residual u, physical dual residual v and their algebraic mismatch u+v. These
should cancel in exact arithmetic; a mismatch above twice the acceptance bound
prevents both norms passing at that iterate.

Separately measures (H T)y-H(Ty), propagates it through the actual physical
current-inertia/Faraday actions, and adds the magnetic metric association gap.
This predicts the inconsistent residual paths; the unexplained fraction is
reported without fitting. Original solver logic, CG/refinement budgets,
physical equations and strict rejection are unchanged. No diagnostic vector is
returned as an accepted state or used to repair a solve. All observations are
recorded outside numeric acceptance components.
