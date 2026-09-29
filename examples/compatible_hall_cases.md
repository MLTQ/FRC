# Compatible Hall verification protocol

Explicit spatial controls retain the original nonlinear Hall snapshot through
nr128, variable density and zero-component cases, independent common quadrature,
and radial refinement at fixed nz256. Coupled cases propagate all three
finite-pressure branches on nr8/16/32 and test opposite axial signs, larger
amplitudes, Hall-free MHD, nonlinear tolerances and fixed-grid time refinement.
The slow branch is evolved for 0.5 radian; intermediate/fast branches for one
period. The different declared durations are never described as equal periods.

Finite-amplitude variable-density, zero-initial-flow states evolve with full
compressible momentum, induction and total energy through nr64. Prescribed
scalar-resistive and stronger-density controls remain explicit. Separate
same-grid time cases use a finer numerical time reference. No case is a true
vacuum interface or proof of FRC formation. Cold electrons remain declared.
