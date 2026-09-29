# Independent Hall / finite-mass verification

Checks both helicities and resistance directly against original electron/ion
momentum (residual scaled by individual forces before cancellation), finite-mass cyclotron asymptotes, pair-mass cancellation, and continuum
cylindrical Maxwell derivatives. Actual FE cylinder solutions converge in
physical-volume field error. A distinct parallel Cartesian limit isolates
second-order time error without claiming to be an axisymmetric cylindrical mode.
Physical quadrature independently checks all three energies and midpoint local
heat; a superposed mean toroidal flow checks total angular momentum. A purely
poloidal perturbation generates toroidal field and flow with Hall and does not
in the MHD control. Invalid inputs/real-field symmetry fail, stationary modes
remain stationary and trial stepping leaves caller fields unchanged.
These tests verify the documented linear reduction only.

Additional checks expose the single-large-step 21.46% phase error despite
roundoff energy conservation, verify the independent full-spatial exponential,
check the continuum local heat integral by direct time quadrature, and confirm
signed-bias reversal exchanges helicities in both physics and actual evolution.
