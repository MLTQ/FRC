# Gauss volume and boundary operator

Positive Gauss-Legendre volume nodes/weights, the exact nodal polynomial
derivative D, and boundary interpolation R give Q=W D and Q+Qᵀ=Rᵀ B R.
The decoupled operator is

Q_N = [[Q−RᵀBR/2, RᵀB/2], [−BR/2, B/2]], B=diag(−1,+1).

It acts on volume plus two boundary traces. The lift is [I,Rᵀ]. Q_N has
zero row sums and only the two signed boundary entries in its symmetric part.
This construction follows the framework of
[Chan](https://arxiv.org/abs/1708.01243); geometry/thermal/centrifugal work is
derived separately for our physical state. There are no volume points at the
axis, so volume divisions by r never require a numerical radius floor.
Boundary values are entropy-variable traces, not additional evolved states.
