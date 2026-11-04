# hybrid3d_circuit_check.py

Independent physical check of the 3D bank-load integral, using two coaxial
circular filaments rather than a second copy of the formation discretization.
Run `python -m examples.hybrid3d_circuit_check`.

- `loop_edges`: analytic Biot–Savart vector potential on Yee edges.
- `reciprocity`: the curl of a plasma-loop potential supplies B. The new
  current-to-coil integral is compared with I_plasma times the exact mutual
  inductance M=μ0√(ab)[(2−k²)K(k²)−2E(k²)]/k. Both solid and annular
  chambers are refined through 20, 10 and 5 mm.
  `project_source=True` includes the production `VacuumCoilBasis` correction
  before evaluating linkage against the same independent analytic mutual
  inductance. This more expensive source solve is optional in the fixture.

The source loop is outside gas, the plasma loop is inside it, and the exact
result is independent of discretization. This checks the new load formula;
it does not establish formation, electron heating or FRC stability. Output is
ignored experiment evidence, not a unit test.
