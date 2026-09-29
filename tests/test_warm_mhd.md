# Moving ideal electron/ion tests

The warm path reproduces the existing ideal fluid state bit for bit at both
orders while evolving electron thermodynamics independently. Constant specific
entropy with nonzero radial and axial compression verifies the analytic
pe_dot=gamma pe ndot/n response and its invariant through actual RK stages.
Species energies sum to the represented conservative thermal energy, and the
reference density scale changes coordinates only. A rejected final species
trial leaves all accepted state, time and diagnostics untouched. Invalid
partitions and as-yet unsupported resistive heating are rejected explicitly.

Independent continuum contact/simple-wave convergence and reference-equation
tests are in `test_warm_reference.py`; algebraic conservation alone is not
evidence of correct thermal transport.
