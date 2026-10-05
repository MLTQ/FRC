# Gauss proposals from stored Newton products

Keeps the actual equations, full rounded-stage and endpoint acceptance,
positivity, strict induction checks and chord/Newton fallbacks of ModalWallGauss.
Previous collocation-polynomial extrapolation is only a starting guess.
A multisecant approximate inverse uses up to 60 centered full-Jacobian products
already evaluated by Newton, rather than small noisy accepted-iterate changes.
The original finite-difference Jacobian scale and GMRES settings are retained.

Reuse requires exactly the previous accepted input state/field and same dt.
All proposal memory is copied locally and commits only after successful
endpoint checks. Failure preserves input and previous continuation history.
The least-squares rank truncation affects only the proposal inverse; no physical
Fourier mode, state or RHS is filtered. All histories and saved actual stages
permit independent equation substitution. No physical damping, mass change,
energy repair or accuracy claim follows from reduced iteration counts.

This explicit variant preserves the source of the two preceding measured,
slower acceleration pilots. Accuracy and conservation still require actual
signed continuum-mode trajectories with time/mesh/amplitude controls.
