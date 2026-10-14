# Continuum reference for enclosed-flux driving

Independent constant-eta linear cylindrical diffusion with physical reduced
electron/ion current inertia. The pulse is Bbar=A sin²(pi t/T) for0<t<T,
zero otherwise. Initial fields/current are zero, density fixed, and motion
and Joule-temperature feedback are second order in A. It is not a full
nonlinear or variable-resistance reference.

From Ampere, Faraday and resistive/inertial Ohm, with psi=r² Bbar/2+u:

    (1-de² Delta*) u_t = D Delta* u - (r²/2) Bbar_dot,
    u(0)=u(R)=0, Delta*=d_rr-(1/r)d_r,
    D=eta/mu0, de²=mu/(mu0 n e²), mu=me mi/(me+mi).

Use q_n R=zeros of J1. Then Bz=Bbar+sum b_n J0(q_n r),
Jtheta=sum q_n b_n J1(q_n r)/mu0, and

    b_n_dot + lambda_n b_n = Bbar_dot/[J0(q_n R)(1+de² q_n²)],
    lambda_n = D q_n²/(1+de² q_n²).

The forcing coefficients follow weighted Bessel orthogonality for u:
(r²/2)=sum[-1/J0(q_n R)] r J1(q_n r)/q_n. See the
[NIST Bessel integral identities](https://dlmf.nist.gov/10.22).
Separation/superposition is the standard transient-diffusion method described
by [Haus and Melcher](https://web.mit.edu/6.013_book/www/chapter10/10.6.html);
the enclosed-flux and finite-inertia formulas above are derived here.

`convolution` evaluates the sinusoidal forcing integral in closed form and
its post-pulse exponential relaxation. For fast convergence of physical
current, extract a(t)(r²-R²/2) from B and -2a(t)r/mu0 from J, with
de² a_dot+D a=Bbar_dot/4. Its Bessel coefficients are4a/(q_n² J0(q_n R)).
Only the rapidly converging remainder is summed. This retains the exact
nonzero continuum wall current; it does not impose the production current
basis's zero wall trace. The associated psi polynomial vanishes at both edges.

Inputs are the example's positive physical constants and supplied radial
points. All reference coefficients and fields are independent of the
production mapping, current operator, transport source and time integrator.
The example compares512/1024-mode references at selected pulse phases to
measure truncation before interpreting production differences. This is one
direct mathematical baseline, not a verification framework.

Three full driven runs are recorded in docs/DRIVEN_FLUX_RESULTS.md. Fixed-eta
Bz improves on8→16 cells, but space-time current error remains28.06%→21.14%
and heating7.87%→4.47% low. The evolving-eta case is a physical control
comparison. The wall-current approximation remains an unresolved solver issue.
