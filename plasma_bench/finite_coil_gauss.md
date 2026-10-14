# Full finite-coil Gauss evolution

FiniteCoordinates removes exactly the inactive regular-axis, rigid-side and
end-normal fluid slots. Every nonaxis poloidal flux coefficient evolves; the
toroidal lift obeys side and both end zero-field traces. Fixed mean density,
sound speed and equal-energy magnetic scales define the nonlinear norm. There
is no periodic endpoint, projection of accepted states or coordinate floor.

FiniteInverse supplies a uniform finite acoustic block and full magnetic
resistive/current-inertial block for correction proposals. Its axial derivative
uses the same physical finite Lobatto interval. Ideal/Hall, changing resistance
and density, pressure/current coupling and exchange remain in the complete
nonlinear residual; their omission from the approximate inverse does not split
the step. This initial inverse is intended for initially resting weak-field runs.

FiniteCoilGauss reuses the complete nonautonomous DrivenModalGauss algorithm,
original bounds, line searches, strict induction solve, positivity and accepted
stage source work. Only active coordinates, correction inverse, source-current
identity and spatial field construction change. No extra physical diffusion.

Seven full24us coil pulses now exercise this coupled path. Independent limits,
spatial/time comparisons, transport flags and remaining exterior-resolution
limits are recorded in docs/FINITE_CYLINDER_RESULTS.md.
