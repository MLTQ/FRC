# Species verification protocol and replay

Defines 14 instantaneous cases (degree 3/5, four meshes, quadrature, density,
uniform-density-contrast and field-grid controls), 8 exact full-transit frozen
current evolutions (mesh, time, current sign and physical ion mass fraction),
and 2 complete spectra/discontinuous-state controls. Both original failed
point-transfer variants remain measurable. Partial and terminal records cannot
be overwritten; failures are saved before the run stops.

The auditor replays every array and non-time metric bitwise on the original
runtime, independently reconstructs physical fields at orders 16/20, and checks
each evolution's final accepted step. The complete source import closure and
code/docs archive accompany raw hashes. Reports of conservation and convergence
remain separate; merely completing the protocol is not a physical validation
gate or a formation result.
Explicit module imports are included in the AST provenance traversal.
