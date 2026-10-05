# Independent coupled implicit checks

The sparse uniform moving-fluid Jacobian is compared with centered differences
of the original nonlinear conservative fluxes at two amplitudes; its error must
decrease quadratically. The coupled physical preconditioner is compared with
the separately derived stationary parent linearization, including all physical
mass/current terms, within the declared strict induction solve accuracy. Its
complex stage inverse is checked by substitution into
the real coupled Gauss and Radau two-stage equations, not their Schur
implementations. Radau's returned state must be its actual second stage.

All three full implicit methods must satisfy their original time equations when
reconstructed from actual returned state/stages. Inputs remain unchanged on
successful wave steps and deliberately exhausted nonlinear contact/invalid
steps. Exact zero velocity defines the separate stationary fixture; a wave's
floating volume mean is not assumed to be exactly zero. The immutable saved
eighth accepted contact state reproduces a former chord-only failure; Newton
must recover under the original nonlinear and induction bounds, with no change
to that input. Complete independent continuum trajectories, time/mesh refinement
and energy are separate frozen studies; these focused tests alone do not prove
long-run accuracy or laboratory formation.
