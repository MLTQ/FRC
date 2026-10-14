# Full moving resistive Bessel decay

Runs all coupled nonlinear fluid, species, Hall, current inertia/stress and
magnetic equations with the new prescribed resistive force and heat source.
Initial density/temperatures are uniform, velocity is zero, bias is zero,
`Bz=A J0(qr)`, `Br=Btheta=0`, `psi=A r J1(qr)/q`, and `J1(qR)=0`.
The exactly constrained initial wall psi is the only endpoint approximation.

Independent linear theory gives

    Bz(t)=A J0(qr) exp(-lambda t), Jtheta=q A J1(qr) exp(-lambda t)/mu0,
    lambda=(eta/mu0) q^2/(1+d_e^2 q^2), d_e^2=mu/(mu0 n e^2),
    released energy=(WB(0)+KJ(0)) [1-exp(-2 lambda t)].

Here mu is the physical reduced species mass. This follows directly from
Ampere, Faraday and resistive Ohm with finite current inertia. Five decay
times are evolved. Lorentz/current stress, heating and pressure-driven motion
are retained; at finite A the nonlinear solution need not equal this linear
target. Amplitude reduction distinguishes that effect from discretization.
This filled periodic wall experiment is not FRC formation or a gas closure.

The scalar eta is 1e-4 ohm meter; all dissipated power is explicitly assigned
to electron heat for this experiment. No electron-ion temperature exchange is
claimed. Geometry is the existing R=0.035m, L=0.24m argon fixture, uniform
n=1e20/m^3, Te=5eV, Ti=1eV with the fixture's gammas 5/3 and1.4.
The gamma_i value is a declared numerical fixture, not monatomic-gas transport.
The constructor offers nr8/16 only, with nr/2 fluid elements and16 axial points.
It uses the original Gauss residual bounds and no artificial diffusion/floor.

Independent root-polynomial/Fourier quadrature20 reconstructs B, J, n, C/n,
velocity and the EOS. Every saved step records field/current error, physical
energy components and physical versus native total-energy drift normalized
by initial analytic magnetic-plus-current energy. Thus the large thermal
background cannot hide a field-energy error. Early and late relative field
errors retain their initial and remaining-signal denominators separately.

Every accepted state and stage, complete configuration, source archive and
solver history is saved. A rejected step remains a terminal failed experiment;
no smaller retry is substituted. The small nr8/16 experiment checks its own
1GiB peak RSS and256MiB disk budget at initialization and each accepted step;
these are discrete observations, not a continuous interrupt monitor.
No unit tests, repeated-protocol framework or automatic replay is supplied.
