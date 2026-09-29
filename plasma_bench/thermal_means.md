# Thermodynamic means

`logarithmic_mean(a,b)=(b-a)/(log b-log a)` uses a symmetric near-equality series
through sixth order to avoid cancellation; arguments must be finite and positive.
Moderate ratios use atanh of the symmetric relative difference, avoiding
subtraction of two large logarithms of tiny SI temperatures. Large ratios use
log1p, with a logarithm-difference fallback only if the ratio itself overflows.
The midpoint is formed from the smaller endpoint plus half the difference,
avoiding needless underflow for equal very small positive arguments.
`entropy_mean(a,b)=a*b/logarithmic_mean(a,b)` avoids the small product directly.
It is the reciprocal logarithmic mean of reciprocals, bounded by its endpoints.

For electron entropy proxy w, entropy_mean times (1/wL-1/wR) equals
log(wR)-log(wL). Thus conservative current flux preserves physical logarithmic
entropy rather than a surrogate quadratic invariant. The temporal logarithmic
mean supplies the exact finite-difference derivative of log U. These are
equation discretizations, not thermal state repairs or fitted means.
