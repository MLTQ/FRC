# Fixed-wavelength radial response and strict refinement

## Purpose and boundaries

The medium full response developed an extra mode concentrated near radial
Bessel15/16. A changing largest numerical eigenvalue cannot establish its
physical meaning. This study fixes the independent radial16 Bessel wavelength
and varies fluid radial elements, field radial cells and axial points separately.
It is an instantaneous consistency study in a filled periodic cylinder with
impermeable fluid wall, fixed axis/wall psi and zero wall Btheta. It is not a
formation calculation or a model of the user's proposed chamber.

The fixture retains R=.035 m, L=.24 m, n0=1e20 m^-3, argon physical ion/electron
masses, Te=5 eV, Ti=1 eV, B0=.005 T, gamma_e=5/3 and gamma_i=1.4. Stationary
n(z), constant pe+pi and the exact ambipolar E are inherited from the documented
contact fixture. These are verification parameters, not recommended lab values.

## Independent parent equations

The independent state is (dn,vr,vtheta,vz,dse,dsi,Br,Btheta,Bz), with J0(qr)
for scalars/vz/Bz and J1(qr) for transverse vectors; J1(qR)=0 at the fixed 16th
root. Let M=mi+me, mu=mi me/M, alpha=me/M, beta=mi/M and a=(-beta,alpha).
The linearized original equations include

- dn_t = -q n vr - partial_z(n vz);
- M n vr_t = q(dp_e+dp_i)+B0 Jtheta;
- M n vtheta_t = -B0 Jr, M n vz_t = -partial_z(dp_e+dp_i);
- ds_s,t = -s_s' (vz+a_s Jz/(e n)), dp_s=gamma_s p_s dn/n+p_s ds_s;
- J=(-partial_z Btheta, partial_z Br+q Bz, q Btheta)/mu0;
- [I+C diag(mu/(e²n)) C/mu0] B_t = -C F + pressure_curl,
  F=-v cross B0+(beta-alpha) J cross B0/(e n), with toroidal
  pressure_curl=sum_s a_s q(p_s' dn-n' dp_s)/(e n²).

The Fourier/Bessel reference imports no production fluid/field/current/time
solver. Its even-grid Nyquist convention is retained and independently tested
by off-grid finite-difference curl and N64/128/256 shape/current/E refinement.
Converged eigenvalue alone is insufficient. Arbitrary eigenvector scalar gain
is aligned only for reference-shape comparison; no evolved state is fitted.

## Actual full RHS comparisons

Sample the same reference into actual conservative fluid and magnetic spaces,
retain the one-time initial wall representation explicitly, and evaluate both
signed nonlinear states. A centered signed difference gives the physical rate.
Reconstruct physical n, velocities, entropy, B, Jt and full E independently at
quadrature16/20. Retain volume/axis/wall/max errors and initial representation.
One-sided exact polynomial traces retain normal-current and tangential-B
errors and jumps at internal faces, including the periodic seam. Strong curl
comparisons exclude surface distributions and use an explicitly stated norm.

No viscosity, resistivity, artificial mass, atmosphere, filtering, positivity
floor or heat/energy repair is added. Remaining numerical dispersion/diffusion
comes from the declared finite spatial approximation and the unchanged strict
linear solve. This rate study has no time-discretization error. It does not
quantify a single effective diffusion coefficient or establish full stability.

## Dual residual refinement

For the admissible magnetic lift T, the original positive induction equation is
T^T[H+Aw^T I(n) Aw]T y = T^T H Fw M^-1 f, Bt=T y,
Aw=-M^-1 Fw^T H. Both preassembled dual action (H T)y and reconstructed H(Ty)
are algebraically equal but can round differently. The saved fine input exposes
opposing operator and physical dual residuals near the unchanged 1e-10 bound.
The original physical-only correction can leave the operator residual failing.

The separate BalancedWallCore uses the mean of the operator and physical dual
corrections when the operator residual fails; otherwise the original correction
is retained. This is a proposal for the next linear iterate. Every original
operator/physical/coefficient/Ohm/Ampere/restored-kernel acceptance check remains,
with the original three inner and five outer refinement budgets. Rejection is
atomic. This changes iterative arithmetic, not the physical operator or its
accepted bounds. Both old passing and failed archives remain reproducible.

Independent polynomial Etheta=x(2-3x²+x⁴), Ez=.37, x=r/R yields analytic
Bzt=(-4+12x²-6x⁴)/R and Jtheta,t=24x(x²-1)/(mu0 R²). With actual nonuniform n,
forcing is the physical weak integral of E-mu Jt/(e²n). Complete reconstructed
fields, -Jt.f=Bt.H.Bt+Jt.I.Jt and atomic exhaustion test the new solver.

## Acceptance and continuation

Require source-frozen full bitwise reruns, all strict residuals, independent
quadrature, amplitude and separate mesh controls. Accuracy remains a separate
gate even when every solve passes. Higher radial/axial resolution, original
local work budgets and full signed trajectory time/residual convergence remain
necessary before moving exact vacuum, transport, coil/end and formation work
can be accepted. Closed flux and axisymmetric equilibrium do not prove
three-dimensional stability or particle confinement.
