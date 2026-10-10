# Kernel-restored linear-inverse wall solve

Separate opt-in wrapper dispatching to SplitInverseWallCore. All original
restored full Ohm/Faraday and kernel reconstruction checks, outer refinement
budgets and diagnostic definitions are copied unchanged from the frozen
balanced solver. Only core electric reconstruction uses linearity of the mass
inverse. Every physical/coefficient bound remains enforced; failure is atomic.
This wrapper is necessary to preserve frozen old source dispatch and evidence.
