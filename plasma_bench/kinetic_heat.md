# Linear kinetic parallel heat transport

An independent uniform-background electron kinetic component, not installed in
the moving formation solver. It evolves a small Fourier temperature perturbation
with full linearized electron-electron Landau collisions and cold, infinitely
massive ion Lorentz scattering. Ions are stationary; density and current are
constrained to zero perturbation by the quasineutral electric field. No sheath,
vacuum edge, changing field direction, finite drift or nonlinear distribution
is represented. The parallel result also applies to a straight unmagnetized
one-dimensional perturbation. It cannot by itself close the actual pulse.

The starting operator is the Landau/Rosenbluth operator in
[Pfefferle et al. (2017)](https://arxiv.org/abs/1701.08037).
Linearization about one Maxwellian gives a symmetric collision bilinear form.
This implementation integrates that form directly, including both test-particle
and field-particle electron terms; it does not use a BGK relaxation surrogate.

Let s=v/vT, vT=sqrt(T/m), T in joules, M(s) a standard 3-D Gaussian,
nu0=n*4*pi*(e²/(4*pi*eps0))²*lnLambda/(m²*vT³), and K=k*vT/nu0.
Time is nu0*t. Backgrounds and logarithm are constant, positive declared inputs.
The dimensionless positive collision matrices are

    Cee_ab = 1/4 <[grad phi_a(s)-grad phi_a(s')]
                      .U(s-s').[grad phi_b(s)-grad phi_b(s')]>_(M M'),
    Cei_ab = Z/2 <grad phi_a(s).U(s).grad phi_b(s)>_M,
    U(g)=(I-ghat ghat)/|g|.

Axisymmetric orthonormal basis: phi_pq=He_p(sz)/sqrt(p!)*L_q(sperp²/2),
p+2q<=degree. Density and parallel velocity are explicit basis elements.
A fixed orthogonal rotation of degree-two elements exposes energy
phi_T=(s²-3)/sqrt(6). Its collision integrand vanishes analytically; density
and electron-self momentum integrands also vanish. These identities are applied
before quadrature, not by correcting a state or fitted collision matrix.

For electron pairs, V=(s+s')/2 is Gaussian with variance1/2, g has variance2.
After the 1/|g| factor the radial measure is exp(-x)dx/sqrt(pi), |g|=2sqrt(x).
Use tensor Gauss-Hermite for V, Gauss-Laguerre for x and Gauss-Legendre for
ghat.z. Axisymmetry permits fixing g azimuth before the isotropic V integral.
Ion radial measure is sqrt(2/pi)*exp(-x)dx, |s|=sqrt(2x). There is no speed
cutoff or small-speed floor. Polynomial degree sets adequate integration orders;
the matrix is accumulated as weighted gradient Gram products.

Multiplication by sz is the exact truncated Hermite recurrence S. On the
density/current-free subspace, da/dtau=(-i*K*S-Cee-Cei)a. Omitting the current
row is the electric-field Lagrange multiplier enforcing zero current, including
its collisional friction contribution. This is an explicit fixture assumption,
not deletion of current inertia or fields from a formation run.

deltaT/T=sqrt(2/3)*a_T; q/(n*T*vT)=sqrt(3/2)*S[T,:]a. Thus
d(deltaT/T)/dtau=-(2/3)*i*K*q. The spatially averaged perturbation energy is
zero for this Fourier mode. Quadratic perturbation free energy decays through
the collision Gram form; no numerical viscosity or damping is added. Truncating
velocity moments can cause recurrence, assessed by increasing degree. Matrix
exponentials are evaluated through eigenmodes; eigenvector conditioning is saved.

The local conductivity follows the odd-moment collision solve q^T C^-1 q,
in units n*vT²/nu0. A separate pure-Lorentz solve is compared with the analytic
Maxwellian speed integral 64*sqrt(2/pi)/Z. Full electron thermalization remains
enabled in the actual temperature-wave evolutions. The example also compares
collisionless evolution with an independently derived Gaussian-characteristic
Volterra solution. No tests or pass/fail framework is required by this module.

The completed saved operators also reproduce the stationary Z=1 leading
nonlocal coefficient from Brodrick et al., Table I: degree28 gives43.5451
against published43.5. The derivation and input-identified one-off calculation
are in outputs/kinetic_heat_published_comparison. This coefficient uses the
stationary non-temperature subspace and each degree's own local conductivity;
it is not a comparison of the time-dependent short-wave modal conductivity.
