# Complete conservative finite-mass trajectory protocol

Eighteen cases evolve forced and unforced nonlinear moving Ar+ states at nr8/16/32,
original-density forced controls, manufactured time refinement, fast compressible
waves through one quarter period, a weaker wave and time controls. The slower
branches have independent rate tests; their very long periods are not silently
called verified by this short fast-wave study.

Every accepted step retains mass, momenta, energy/source work, species extrema,
entropy diagnostics and solve work. Initial/final full states are archived.
Failure preserves the exact last accepted input, attempted dt and exception;
unchanged fields and budgets are checked and no hidden retry occurs. A partial
trajectory remains partial. Error metrics are computed at its actual time.

Source/configuration/parameter/field hashes freeze each run. A read-only audit
must verify histories, independent physical energy and final errors before a
collection is declared complete. Global energy is an algebraic consequence of
shared face fluxes, not evidence that species temperatures or slow modes are
accurate. The prototype explicit RK4 cap requires measured time sensitivity.
