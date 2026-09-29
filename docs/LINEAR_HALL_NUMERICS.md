# Axisymmetric linear Hall / finite-electron-mass reference

This component evolves both poloidal and toroidal magnetic field and bulk
velocity in a fixed, uniformly filled cylinder. It is a controlled linear,
incompressible reference about B0 ez, n=constant, vi=ve=0. It is not a nonlinear
formation solver or the free plasma/vacuum boundary. The original low-pressure
argon target, heated-edge problem and Stage 2 guard remain unchanged.

## Parent equations and work

For singly charged ions and electrons, let M=mi+me, rho=n M,
v=(mi vi+me ve)/M, J=en(vi-ve), mu=mi me/M. Linearizing species momentum gives

    mi n d_t vi = en(E+vi cross B0) - grad pi - en eta J
    me n d_t ve = -en(E+ve cross B0) - grad pe + en eta J
    curl b = mu0 J, div b=0, d_t b=-curl E, div vi=div ve=0.

There is no displacement current. Uniform density makes the species pressure
combination in Ohm's law a pure gradient. Summing and eliminating species
velocities, then taking the solenoidal projection P, gives

    rho d_t v = P(J cross B0) = B0/mu0 d_z b
    E + v cross B0 = alpha J cross B0 + eta J + beta d_t J + gradient
    alpha=(mi-me)/(M en), beta=mu/(n e²), d*²=beta/mu0
    (1+d*² curl curl) d_t b = B0 d_z v - alpha B0 d_z J - eta curl J.

The pressure constraint cancels the magnetic-pressure gradient; this is an
incompressible reduction, not a cold compressible oblique wave. In the exact
parallel transverse limit the pressure gradients vanish. Current drift energy
is beta J²/2. Indeed the sum of species kinetic energies equals
rho v²/2 + beta J²/2. Integrating the work gives

    d_t integral [rho v²/2 + b²/(2 mu0) + beta J²/2] = - integral eta J².

The Hall force has zero J-dot work; magnetic/bulk transfer cancels. Heat is
accounted as positive dissipation, not fed into a thermal equation in this
linearization. No unexplained energy repair is used. Physical electron mass
is retained exactly in `two_fluid`. `hall` sets me=0 consistently in all
coefficients; `mhd` additionally sets alpha=0. These comparisons have explicit
names. Changing the ion mass is never a hidden electron-mass modification.

These equations are derived here from the species system. Primary context:
[Andres et al., species equations and kinetic-energy invariant](https://arxiv.org/pdf/1410.6518),
equations 1–6, 10–12 and 26. Their nonlinear calculations and transport choices
are not validation of this component. Nonlinear advection, drift stresses,
species thermal transport, neutrals, finite-Larmor-radius and kinetic response
remain absent. Finite inertia alone does not justify a fluid closure at low n.

## Cylinder, axis and boundary

Use b=curl(a e_theta)+t e_theta and v=curl(u e_theta)+w e_theta.
The poloidal flux is psi=r a (the uniform background adds B0 r²/2).
Regular a,t,u,w vanish on axis. At r=R set a=t=u=w=0, with periodic z of length L.
These homogeneous reference traces set br=Jr=vr=0 and btheta=vtheta=0.
They give zero perturbation energy flux and no integrated axial torque.
They are NOT advertised as a perfectly conducting wall: Ez can be nonzero
with resistance. They impose a restricted zero-flux perturbation sector and
are not a quartz/vacuum interface. There is no exterior computational gas.

Radial piecewise-linear nodal elements with zero endpoint values use

    M_ij = integral_0^R r phi_i phi_j dr
    S_ij = integral_0^R r (phi_i'+phi_i/r)(phi_j'+phi_j/r) dr.

Positive Gaussian quadrature resolves the 1/r term away from the axis;
on the axis element the sole active shape is proportional to r and regular.
Generalized diagonalization S Q=M Q Lambda, Q^T M Q=I diagonalizes this
numerically assembled operator; no analytic Bessel root enters production.
Axial differentiation uses the Fourier symbol ikz on an odd periodic grid
(no ambiguous Nyquist derivative). Every radial degree of freedom is retained.
Input data must resolve its axial spectrum; sampling cannot detect aliasing.

For each radial/Fourier mode K²=Lambda+kz², normalize poloidal coefficients
by K, velocity by sqrt(mu0 rho), and write x=(Bp,Bt,Up,Ut). Curl swaps the
poloidal/toroidal pair with K. The time mass is diag(A,A,1,1), A=1+d*² K².
Ideal velocity/magnetic couplings are ikz vA; Hall magnetic off-diagonals
are -i alpha B0 kz K/mu0; resistive diagonals are -eta K²/mu0.
The ideal operator is skew Hermitian in this work metric. Implicit midpoint
therefore preserves quadratic energy to roundoff for eta=0 and loses exactly
its separately integrated midpoint eta J² for eta>0. This is a property of
the numerical time scheme, not evidence of phase accuracy. No artificial
resistivity, smoothing, limiter or energy redistribution is included.

The recovered toroidal current is the weak FE curl, not an unreported strong
second derivative of a piecewise-linear potential. Poloidal fields/current
use actual first derivatives, with the regular axial limit Bz=2 a'(0).
Physical-volume quadrature independently checks magnetic, velocity and weak
current energy and local heat. The angular momentum density is rho r vtheta.

## Independent continuum waves

Let q=j_(1,m)/R be a positive zero of J1, k=2 pi ell/L, K=sqrt(q²+k²),
and h=+1 or -1. Set

    a=C J1(qr) exp(i kz-i omega t), t=h K a,
    br=-ik a, bz=C q J0(qr) exp(i kz-i omega t), curl b=h K b,
    v=-k B0/(mu0 rho omega) b.

Both helicities are regular at the axis and satisfy the reference boundary.
The dispersion, derived by substitution independently of the matrices, is

    A omega² + (i D K²-H) omega - (k vA)²=0,
    A=1+d*²K², D=eta/mu0, H=h alpha B0 k K/mu0.

Positive-frequency roots for k,B0>0 distinguish the whistler/electron branch
(h=+) and ion branch (h=-). Cylindrical polarization is generally elliptical.
[Krishan and Varghese's author manuscript](https://www.researchgate.net/publication/225944123_Cylindrical_Hall_-_MHD_Waves_A_Nonlinear_Solution),
equations 23–26, provides the cylindrical curl eigenfunctions and the massless
Hall limit in a different, partially ionized rotating setting. Our coefficients
and finite-mass quadratic are independently derived above. No neutral closure
from that setting is borrowed.

For q=0 the separate Cartesian transverse circular-wave limit satisfies the
original species momenta directly with E=i omega b/(h k). At increasing k,
the exact finite-mass branches tend to eB0/me and eB0/mi; massless Hall has an
unbounded whistler frequency. Independent checks include both helicities,
damping, species residuals, mesh/time refinement and energy. An accurate
wave test does not establish nonlinear formation or edge validity.
