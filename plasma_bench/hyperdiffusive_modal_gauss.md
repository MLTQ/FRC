# Implicit proposals including numerical diffusion

Uses the original GuardedModalGauss step, full nonlinear stages, acceptance
bounds and continuation. Only the uniform-background proposal inverse changes.
For a Fourier mode with admissible magnetic lift T, let A=-M^-1 F* H T,
Hr=T* H T, D=A* M A and G=Hr+mu/(e^2 n0) D. The magnetic positive Laplacian
is L=mu0 Hr^-1 D. Add -eta6 G^-1 D L^2 to the original magnetic block.
This is the linearization of HyperdiffusiveWallPlasma at uniform density and
zero current. Numerical heating is quadratic, so adds no such linear block.
No full-grid dense matrices are built. The repeated small-matrix construction
keeps frozen ModalWallBlocks sources intact.

Gauss remains fourth order and A-stable, not L-stable: very stiff decaying
modes can be weakly damped per step. Stability/accuracy cannot be inferred from
this proposal or its solve residual. Assess actual complete evolution against
independent fields and report the numerical dissipation explicitly.
