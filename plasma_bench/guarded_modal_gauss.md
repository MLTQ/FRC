# Bounded proposal work for full Gauss equations

Retains RecycledModalGauss physical equations, previous-step starting guess,
Newton-product memory and all unchanged stage/endpoint/induction acceptance
bounds. Approximate secant/chord proposals get one full-step trial and must
halve the residual (or already pass its original bound). After four nonlinear
iterations only a useful secant proposal or full Newton is attempted. Newton
keeps its original residual-decreasing 13-trial line search.

This prevents slow proposal stagnation from consuming the nonlinear budget.
It tightens only the usefulness criterion for an approximate proposal; it never
relaxes acceptance, retries time, filters fields or repairs energy. Failure is
atomic, and continuation memory commits only at a passing endpoint. Histories
expose every accepted method and rank. Original slower pilot sources remain
frozen. Actual time/mesh/amplitude and physical audits are still required.
