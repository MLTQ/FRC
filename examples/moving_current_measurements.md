# Manufactured moving-current rate protocol

Nineteen cases cover physical Ar+/electron masses at n0=1e18 on nr=8..128,
target n0=1e20 and consistent massless controls on nr=8..64, a physical pair-mass
control, two extra compression phases, higher quadrature and tighter residuals.
Main cases request 1e-10 equation tolerance; separate nr16/64 controls request
1e-12/2e-12. The production rate solver's standalone default remains 2e-12.

Every case retains exact initial fields/kinematics/density rate and calculated
Bdot/Jdot/full E. A linear rejection retains unchanged inputs and the actual
failure instead of a completed label. Physical norms and independent changing-
density energy derivatives retain explicit linear-refinement counts and all
residual components. The original failed source/data are separately preserved.
One nr32 case separately removes
each density, vorticity and bulk term to expose its effect; those altered rates
are controls, not alternate physical models. The reference includes explicit
parent-species manufactured forcing, not fitted discrete forcing.

Source is frozen across each measurement, NPZ hashes are saved, and cases can
run independently via --group. A separate read-only report must audit every
expected case before declaring the collection complete. No trajectories or
formation conclusions can be inferred from this instantaneous rate study.
