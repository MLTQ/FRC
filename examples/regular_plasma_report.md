# Coupled verification report

Reads terminal audited raw records and plots measured contact/fast-wave rates,
full-period signal errors, magnetic half-step differences, slower-branch errors
and all stationary spectrum eigenvalues. It includes the difficult slow branch
and the remaining strict-residual/contact/interface/vacuum gates explicitly.

Magnetic temporal differences use the fixed physical magnetic Hessian on
successive final solutions, normalized by the resolved perturbation. The plot
does not mistake the background field for signal or spatial error for temporal
error. Saves machine-readable summary, PNG and PDF. This is presentation only;
it does not modify frozen source/results or change any accepted state.
