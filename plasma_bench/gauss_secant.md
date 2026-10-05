# Full-residual multisecant proposal

For a fixed approximate inverse P and accepted nonlinear pairs S=delta x,
Y=delta r, G=P Y, propose dx=-P r-(S-G) least_squares(G,P r). This satisfies
B Y=S in independent secant directions. Each pair is normalized by ||P delta r||;
relative singular values below 1e-10 are omitted only in the proposal fit.
No physical field, rate or accepted state is filtered. Full-equation residual
checks and line search determine acceptance in accelerated_modal_gauss.

Stores at most 16 pairs by default, local to one time step. A failed time step
cannot leave cross-step history. The fixed inverse must stay unchanged during
its lifetime. Zero/nonfinite images are ignored; independent dimensions and all
physical modes remain available through the base inverse and Newton fallback.
