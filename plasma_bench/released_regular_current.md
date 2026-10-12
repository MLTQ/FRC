# Release expired electromagnetic workspace

`ReleasedRegularCurrent` inherits the original actual-density current model.
Its electromagnetic pair retains every arithmetic expression and call order,
releasing sampled current/velocity/magnetic arrays after ideal/Hall evaluation
and the physical force array after its work-adjoint fluid rate is computed.

Only local references expire: no caller array, state, matrix or rate is changed.
Original density forcing and explicit incomplete density exchange are retained,
along with all actual physical masses and operator definitions. No new numerical
diffusion, atmosphere, floor, filtering, energy repair or acceptance rule.

Require source AST identity after stripping added del statements, followed by
complete original/new signed physical workflows, exact input/rate/audit
comparisons and a fresh full replay. Measure own-process peak storage instead
of asserting a reduction from theoretical array counts. The original model
remains unchanged; this subclass is selected only in the opt-in workflow.
