# Accelerated full Gauss nonlinear solve

Opt-in subclass of ModalWallGauss. Solves exactly the same two-stage Gauss
collocation equations, with identical original rounded-stage, endpoint,
positivity and induction residual checks. Constructor, physical coordinate
scales and uniform modal inverse are inherited. No time retry or tolerance
relaxation is introduced.

Before the ordinary chord/Newton fallbacks, a GaussSecant proposal reuses
accepted nonlinear iterates. Every proposal has the same full-equation
residual-decreasing line search. Histories retain accepted method, secant rank,
backtracking, Newton work and every original check. Memory is local to a step;
failed solves preserve the supplied state/field and carry no history forward.

The explicit step implementation is copied from the frozen original Gauss
algorithm to keep previous source closures immutable. The only algorithmic
addition is the proposal and its diagnostic counters. Tests independently
substitute saved actual stages and endpoints into the original equations.
This changes solver work, not temporal order, physical diffusion or physics.
