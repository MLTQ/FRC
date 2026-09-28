# Independent coupling fixtures

`aligned_meshes` provides a unit SI fluid cylinder and a larger aligned magnetic
domain. `compact_flux` is a C3 polynomial packet that vanishes outside its support.
`compact_force` differentiates its continuum polynomial analytically to obtain
Jtheta and both components of J cross B, including an optional uniform bias.
`cell_average` uses independent tensor Gauss quadrature and cylindrical weights.
These exact force averages do not reuse the numerical stiffness or interpolation.

`passive_advection` initializes a very weak flux packet in uniform axial flow.
The field force is O(amplitude²); the analytic packet translates at the prescribed
speed until boundary-driven gas waves can affect its support. Closed ends are
still enforced and global gas energy includes the resulting end disturbances.
`gaussian_diffusion` solves unbounded psi_t=D Delta*psi exactly: a regular r²
Gaussian broadens with width²+4Dt and amplitude factor (width²/spread)^(5/2).
`passive_diffusion` initializes its harmonic-vacuum continuation. Comparisons are
restricted to an interior window before its exponentially small boundary tail
matters; the finite chamber is not claimed to be the exact unbounded solution.

`coil_pulse` is a fixed-chamber integration fixture with two analytic filaments,
uniform bias, scalar resistivity and a smooth unipolar current pulse. Amplitude
is normalized to 0.3 times the bias at the chamber center in vacuum. Its unit
geometry and plasma constants exercise coupling, not experimental argon transport.
The pulse time step is capped at 0.002 s (at least 20 steps per rise/fall), in
addition to CFL/diffusion constraints; temporal refinement is reported separately.

All dimensions, times and fields are SI verification fixtures, not suggested
experimental parameters. Cell force quadrature evaluates strictly positive r;
`compact_force` itself is not an axis-point evaluator. These smooth fixtures do
not establish shock fidelity or the validity of a free plasma boundary.
