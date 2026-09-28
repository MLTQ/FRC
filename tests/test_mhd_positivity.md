# Admissible reconstruction verification

The density test supplies positive cells whose combined 2-D slopes previously
made a corner negative. The thermal test uses the independently derived curved
equilibrium and checks every flux/UCT trace without modifying the input state.
Higher-order independent quadrature verifies that the radial force source uses
the same modified total pressure as the reconstruction.

Resolved wave evolution must be bit-for-bit unchanged with limiting enabled.
The formerly rejected p0=0.1 Pa equilibrium must finish on 12/24/48-cell grids
with conserved mass/energy, div B and decreasing drift; limiter activity must
decrease to zero on the finest tested grid. Rescaling B and pressure checks
that no absolute SI pressure floor has been introduced. Pure RHS queries must
not change accepted-stage counts. The explicit disable switch retains reference
failures and does not imply that disabling the limiter is physically preferable.
An independently solved magnetized double rarefaction checks actual evolving
low-density/low-pressure states, conservation and unchanged parallel B. Its L1
errors must decrease under refinement. The generated report separately retains
the more difficult central pressure error; passing this test is not a claim
that all rarefaction details are resolved on the tested grids.
