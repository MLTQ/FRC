# Independent Gaussian-characteristic slab

## Purpose

Analytic free-streaming comparison, separate from numerical advection and its
moment routines. Integrate f=phi(s)*1(|x-vT*s*t|<a) over each actual spatial
cell and over the declared finite velocity interval [-sigma,sigma]. This uses
no numerical trajectory, reconstructed flux, moment stencil or fitted field.

## Method and contract

`cell_moments` returns raw velocity moments k=0..3, averaged over x. At fixed s,
the occupied fraction is length([xl,xr] intersect [-a+vT*s*t,a+vT*s*t])/dx.
Its exact breakpoints are (xl±a)/(vT*t),(xr±a)/(vT*t) and velocity endpoints.
Between them the fraction is linear. Integrals of s^k phi(s) follow
I0=Phi(u)-Phi(l), I1=phi(l)-phi(u),
Ik=(k-1)I(k-2)+l^(k-1)phi(l)-u^(k-1)phi(u).
Tail CDF differences use the negative tail to avoid subtracting numbers near1.
The resulting cell averages need no spatial quadrature or sampled edge masking.

| Dependent | Expects | Breaking change |
| --- | --- | --- |
| `examples/vlasov_slab.py` | dimensionless cell-average M0..M3 | normalized/truncated moments or changed support |

Time must be positive. This is a mathematical free-transport limit, not an exact
solution of self-consistent Vlasov-Poisson or the isotropic fluid model. The
velocity tail is deliberately identical to the declared numerical domain;
its omitted probability/energy must still be reported separately.
