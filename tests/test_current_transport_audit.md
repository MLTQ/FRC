# Independent transport-audit checks

Verifies root-polynomial values, first/second derivatives and cell-edge average
reconstruction against analytic polynomials, including periodic translation.
Checks the declared variable-density physical velocity derivative using centered
finite differences and Richardson cancellation, separately from production
adjoints or energy assertions. These tests support the audit machinery; the
sixteen measured cases establish spatial convergence and expose conservation
limits independently.
