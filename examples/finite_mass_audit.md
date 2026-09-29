# Frozen conservative-rate and actual-trajectory audit

Validates the baseline source archive and every old source against current code,
allowing only the explicitly expanded provenance-driver roots to differ. Each
record's configuration, parameters, NPZ hash, initial arrays and source identity
are checked. Rate audits replay conservative flux divergence, full-E Faraday,
independent exact species targets and order-8/12 norms. Large errors remain data.

Evolution audits inspect every accepted budget/positivity/time entry, endpoint
integrals and physical component energies/entropy. Independent analytic parent
total-energy change checks manufactured source work, with explicit RK quadrature
error retained. Initial and final representation errors and tiny perturbation
normalization remain distinct. Rejected trajectories stay partial.
Unforced axial/angular momentum histories are checked against the initial
integral of their absolute local density, so cancellation of net wave momentum
does not produce an invalid relative denominator. A first-step failure can have
only the unchanged initial history row; it is never accepted as a completion.

This audit is read-only. It does not independently rerun every accepted state;
the complete source, histories and exact final fields identify those trajectories.
Analytic oracles establish numerical accuracy, while archive/history checks
establish integrity; neither alone proves formation or physical closure validity.
