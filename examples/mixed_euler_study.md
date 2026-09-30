# Mixed radial and axial evolution protocol

The identical smooth mixed initial state is evolved either with its analytic
continuous manufactured source or with no forcing. The manufactured target
is exact by construction; unforced evolution has no analytic final target.
Mesh refinements2/4/8/16 elements and16/32/64/128 axial points accompany
fixed-grid step halving/quartering and separate radial/axial resolution controls.
Duration is one initial radial sound-crossing unless explicitly recorded.

All physical states use the frozen cylindrical solver. The manufactured RK
wrapper integrates independent source work and conserved source totals only
as diagnostic ledgers. Actual energy-work and mass/angular/axial/entropy
source remainders are saved; the ledgers never correct the state.

Final independent physical comparisons use a common32-element,24-Gauss-node,
256-point axial grid, with a fixed central region r<R/8. They retain both
species, all velocities, maximum errors and actual continuous energy/integrals.
Higher quadrature and all saved snapshots are checked by the separate audit.

Atomic accepted checkpoints are written every roughly10 percent. Terminal
failures preserve accepted state and rejected step. The source is frozen before
running; terminal records cannot be overwritten. Independent cases can use
local processes; this is numerical parallel execution, not agent delegation.
No resume is implemented; never restart a live or terminal failed case because
an observation timed out. These are fluid verification fixtures, not formation.
