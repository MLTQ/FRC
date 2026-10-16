# Velocity resolution on the failed moving reacting states

This bounded numerical study assesses a concrete proposed solver change:
keep the exact fine density/support partition while using fewer mechanical
velocity degrees of freedom. It uses the actual first accepted states of
`reactive_material_argon_collapse`, including all occupied ground-Ar charge
states0..18 and electron index19, physical masses, and unchanged0-frame moments.
It does not alter a source dataset, advance a trajectory or certify a solver.

`restrict_velocity` selects a nested cosine-distributed subset of density faces
with the original requested number of velocity cells. All physical support ends
remain. Analytic polynomial integrals assemble M_ij=integral(rho phi_i phi_j dx)
and b_i=integral(rho v phi_i dx). M v_new=b is the weighted L2 projection. Nodal
partition of unity preserves total vector momentum. Fine-cell thermal recovery
uses U_new=U_old+K_old-K_new with the actual exactly integrated P1 kinetic energy.
Any negative recovered energy rejects the candidate; none is clipped or repaired.
Global kinetic loss can be positive while individual cold cells are inadmissible.

Density, charge and geometry are unchanged, so the represented Gauss field is
unchanged. This is the attraction of separating quadrature cuts from velocity
degrees of freedom, but does not establish admissible local kinetic/thermal
resolution or convergence. Source-generated velocity structure still matters.

The adiabatic acoustic stiffness restricted to the same velocity space is
K_ij=integral(gamma p dphi_i/dx dphi_j/dx dx). Its largest generalized eigenvalue
K u=omega² M u supplies the recorded acoustic angular frequency. The pressure is
the original state pressure, even when projected thermal recovery fails; the
spectrum describes a candidate space only. Electrostatic and other stiffness
are excluded, so this is explicitly not a full-system timestep bound.

`run` retains every input hash, source archive, projected state and matrix. Use
`python -m examples.material_velocity_resolution SOURCE_OUTPUT FRESH_OUTPUT`.
Only actual saved failure states are examined; no unit tests or test framework.
