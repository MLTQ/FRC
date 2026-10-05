# Independent stationary two-species continuum linearization

This operator does not import production fluid, projection, current, induction
or time solvers. It linearizes their original continuum parent equations about
n(z), pe(z)+pi(z)=constant, v=J=0, B=B0 ez and the exact ambipolar electric
field. Scalars/vz/Bz use J0(qr); transverse vectors use J1(qr), with J1(qR)=0.
Axial derivatives use an independently constructed Fourier collocation matrix.
Even-grid real Nyquist differentiation is zero; refine the axial mesh.

State is (dn,vr,vtheta,vz,dse,dsi,Br,Btheta,Bz). Let M=mi+me,
alpha=me/M, beta=mi/M, a=(-beta,alpha), mu=mi me/M and
dp_s=gamma_s p_s dn/n+p_s ds_s. Linear equations are
dn_t=-q n vr-partial_z(n vz), M n vr_t=q dp+B0 Jtheta,
M n vtheta_t=-B0 Jr, M n vz_t=-partial_z(dp), and
ds_s,t=-s_s' (vz+a_s Jz/(e n)), with s_s'=p_s'/p_s-gamma_s n'/n.
Ampere gives (Jr,Jtheta,Jz)=(-Btheta',Br'+q Bz,q Btheta)/mu0.

Let C be that cylindrical curl. Induction solves
[I+C diag(mu/(e² n)) C/mu0] B_t=-C F+pressure_curl,
F=-v cross B0+(beta-alpha) J cross B0/(e n).
Only the theta pressure-curl component is nonzero:
sum_s a_s q (p_s' dn-n' dp_s)/(e n²). This differentiates the full
pressure force analytically, including the density variation multiplying
the background pressure gradient. No compatible projection or numerical
curl-kernel subtraction from the production solver is used.

The three original uniform finite-mass Bessel/Fourier wave branches provide
independent sign, mass and pressure-gradient controls. A variable-profile
spectrum remains a continuum *discretization* requiring axial/radial mode
refinement, divergence checks and independent perturbation evolution. A
positive mode must not be assigned to physical plasma without those controls.
