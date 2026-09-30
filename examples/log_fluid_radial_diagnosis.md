# Radial acoustic stability diagnosis

Linearizes the conservative continuity and pressure/hoop maps about uniform,
stationary gas at constant species entropy. Dimensionless density and radial
velocity give A=c[[0,-D],[G,0]]. Independent central differences of the full
nonlinear fluid RHS verify that this block is its actual acoustic limit.

The true physical density and radial-momentum metrics define the diagnostic
energy-adjoint gradient Hm^-1 D^T Hn. Its spectrum is compared, without adopting
it in the nonlinear solver or claiming full conservative momentum/energy.
The fastest unstable eigenvector is evolved by an independent matrix exponential
and RK4 at three steps to separate spatial growth from time error. Degree3/5,
nr8/16/32 and higher quadrature retain matrices, complex eigenvectors and all
measurements. No damping or fitted physical parameter is introduced.
