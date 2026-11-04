# iec_view_fields.py

Read-only geometry for native electrostatic field lines and an optional magnetic
diagnostic. It never changes the solver, particle state or physical field.

- `electric_lines`: midpoint arclength integration of ±E/|E| in the saved
  piecewise-constant FEM field. Native mesh face tracing stops at absorbing
  steel. Trace spacing is a display integration choice, not new field resolution.
  Arrow order follows +E. Seeds are independent full-sphere display locations;
  no rotated copies or assumed plasma symmetry.
- `biot_savart`: μ0/(4π) sum (J dV)×R/(R²+ε²)^1.5. Nodal current moments
  include every saved physical macroparticle. Softening is explicitly the bulk
  mesh target; this is a regularized plasma-only contribution, not total B.
- `magnetic_lines`: 17³ Cartesian sampling and midpoint trace, inside the
  chamber. The field ignores feed/return/electrode/displacement currents and
  steel permeability. Open transient plasma currents do not define a complete
  Ampère/Maxwell solution. No Lorentz force or induction is added to the run.

Analytic comparisons use uniform-potential-gradient lines and the on-axis field
of a closed circular current loop. A magnetic picture cannot establish trapping
or justify the electrostatic approximation on its own.
