# Variable-density component protocol

Runs degree 3/5 at nr8/16/32/64 with nz=2nr, plus nr32 quadrature8/12,
physical density1e16/1e20, reduced density contrast, and reversed continuity
motion controls: **14 instantaneous cases**. The default is n0=1e18 and physical
Ar+/electron masses. Actual integrated prescribed particle flux supplies n_t;
an independent volume integral of the analytic divergence checks it.

Every case solves the original moving-density induction with explicitly declared
manufactured source, saves all input/output cochains and weak forces, and checks
local and global energy derivatives. Independent primitive/root-product norms at
quadrature12/16 report n/n_t/B/J/B_t/J_t/E, pressure, energy and axis/max errors.
Both gammas retain separate graph work/log-entropy and nodal entropy-rate errors.
The erroneous old-stream work is an explicit omission control, not adopted.

Fresh output directory only. Freeze source before cases; preserve exceptions with
their traceback rather than silently retrying or relabeling initialization.
Any failed case or changed source makes the process fail. This is a bounded
manufactured component study, not a time trajectory or complete two-fluid result.
