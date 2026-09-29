# Bias precision diagnosis

Uses the original rounded nr32 slow case, its exact stored velocity and the
unchanged numerical field operators. Compares its original current/rate with
subtracting the known harmonic after rounding and independently retaining the
analytic initial perturbation before any addition. The latter is a representation
control, not a production exact-current/velocity substitution. All three input,
current and rate arrays plus the original input/source hashes are retained.
Also measures the false current obtained by differentiating only the uniform
harmonic through the rounded matrix action. No earlier case is overwritten.
