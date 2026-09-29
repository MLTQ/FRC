# Nonlinear compressible axisymmetric toroidal MHD

`ToroidalMHD` is a separate six-variable extension of the bounded resistive MHD
core, with all three velocity and magnetic components. It retains compressible
pressure evolution and nonlinear advection, magnetic tension, rotation and hoop
stress. The old four-variable solver and its API are unchanged. This component
is massless-electron MHD; nonlinear Hall/finite-mass transport is still absent.
It is not the plasma/vacuum formation interface or a claim of model validity
for the original low-pressure argon target.

## Equations and annular state

Axisymmetric SI variables are U=(rho,rho vr,rho vz,E,L,g), where
L=r rho vtheta and g=Btheta/r. Corner psi still gives divergence-free Br/Bz.
For poloidal normal direction n and total pressure pt=p+B²/(2mu0), fluxes are

    F_rho = rho vn
    F_mr = rho vr vn - Br Bn/mu0 + pt delta_nr
    F_mz = rho vz vn - Bz Bn/mu0 + pt delta_nz
    F_E = (E+pt) vn - (v dot B) Bn/mu0 + (E_res cross B)_n/mu0
    F_L = r (rho vtheta vn - Btheta Bn/mu0)
    F_g = g vn - (vtheta/r) Bn
    d_t U = -div_p F + (0, [pt+rho vtheta²-Btheta²/mu0]/r, 0, 0, 0, 0)
    psi_t = -r Etheta, Etheta=vr Bz-vz Br+eta Jtheta.

The toroidal resistive contribution is dg/dt=(1/r)d_r Ez-d_z(Er/r).
This follows directly from Faraday, with Jr=-d_z Btheta/mu0 and
Jz=d_r(r Btheta)/(mu0 r). Applying ordinary cylindrical divergence to Btheta
itself would produce the wrong geometric term. L is conservative: no discretized
azimuthal torque source is needed. Total energy includes all three components;
eta J² is diagnosed but never added again to total energy.

Primary context for the cylindrical conservation laws is
[Skinner and Ostriker (2010)](https://arxiv.org/pdf/1004.2487), equations 9–17.
The implementation and regular moment convention here are independently derived;
it is not Athena or a reproduction of that algorithm. The need for actual
curvilinear cell moments is also discussed by
[Mignone (2014)](https://arxiv.org/pdf/1404.0537).

Cell variables are annular volume averages. With s=<r²>=(r_left²+r_right²)/2,
reconstruction assumes local regular vtheta=r Omega and Btheta=r g. Thus

    Omega = L/(rho s)
    E_rotation = rho s Omega²/2
    E_toroidal = s g²/(2mu0).

This is exact for uniform-density rigid rotation and linear Btheta, including
the axis. It is a second-order consistent subcell model, not exact recovery of
arbitrary covariance from cell averages. Exact continuum transport fixtures
initialize true averaged E and L separately: the small resulting pressure
recovery error is retained and measured, not erased by initial energy adjustment.

## Numerical flux, source and boundary

First order uses constant poloidal states; second order uses established odd
vr/Br and even variables in x=r²/2. Omega and g are even. Reconstructed total
pressure includes s g²/(2mu0). Radial traces use actual face radius; axial
quadratic toroidal terms use sqrt(s). The same reconstructed/blended total
pressure enters face flux and radial source. Five-point Gauss integrates the
radial source exactly for these polynomial reconstructions, including degree-
eight rho r² Omega². Density and pressure trace limiters are explicit and
reported; no cell average, density floor or energy reset is introduced.

Local Rusanov jumps dissipate the actual conservative variables with the full
fast magnetosonic speed. UCT also uses the full speed and retains directional
magnetic jumps. These terms are numerical viscosity/diffusion and are not
called physical resistivity. SSPRK2 is subject to a full fast-wave/explicit
resistive bound. Both stages must pass before time, state, heat and limiter
accounting commit. Trace positivity does not prove full-update shock robustness.

The filled cylinder has periodic z and an impermeable, slip, perfectly conducting
outer wall: vr=Br=0, Etheta=Ez=0 there. The axis is a coordinate origin,
not a second conductor; Etheta=0 but regular Ez can be finite. There is no
exterior atmosphere and no vacuum coupling in this bounded reference. Angular
momentum, mass and total energy telescope under these boundaries. Toroidal flux
integral Btheta dr dz can change through Ez at the axis under resistance.

## Toroidal resistive work pairing

Let m_i=integral_cell r dr, s_i=<r²>, D=eta/mu0 and g_i the cell value.
Number radial electric nodes j=0..N (axis to wall) with Ez_N=0. Let

    b_0=s_0 g_0, b_j=s_j g_j-s_(j-1) g_(j-1) for 1<=j<N
    A_j=(j+1/2) dr²/8 for edges 0<=j<N
    W_0=dr²/4, W_j=j dr² for 1<=j<N
    P_jj=W_j-A_(j-1)-A_j, P_(j,j+1)=A_j (A_-1=0)
    P Ez = D b.

The matrix's last diagonal retains A_(N-1), although its wall neighbor is the
prescribed zero field. This compatible current mass satisfies

    E^T P E = sum_j [W_j-2(A_(j-1)+A_j)] E_j²
                + sum_edges A_j (E_j+E_(j+1))² > 0.

All displayed weights are positive. With (Er/r)=-D d_z g, update

    dg_i/dt = (Ez_(i+1)-Ez_i)/m_i - d_z(Er/r).

Summation by parts with EB=pi dz/mu0 sum_i,j m_i s_i g_i,j² gives

    -d_t EB = 2pi dz/eta sum_z Ez^T P Ez
              + sum_cells volume s (Er/r)²/eta >= 0.

The tridiagonal solve is a current reconstruction, not an implicit time step;
SSPRK2 remains explicit and time restricted. Its temporal energy-partition error
is measured separately. The positive square decomposition supplies heat while
an independent magnetic-energy derivative verifies the work identity. Radial
and axial resistive Poynting fluxes enter total energy only once.

The axis current equation is included, rather than imposing Ez=0 at an imaginary
inner conductor. For g=a+b r² with exact cell averages and a=-2bR², the reconstructed
Ez=D(2a+4br²) satisfies the conducting trace and dg/dt=8Db exactly at every cell.
The off-diagonal current mass is essential: diagonal lumping of this metric
would give 10Db in the first cell. Numerical current quadrature is consistent
and positive; it is not a collection of exact geometric dual-cell volumes.
Independent Bessel refinement tests the actual spatial accuracy beyond this
polynomial identity.

## Independent references

Rigid rotation with linear pinch has rho=1, vtheta=Omega r,
Btheta=sqrt(mu0)c r, uniform Bz, p=p0+(Omega²/2-c²)r². Unbalanced constant p
has known radial acceleration (Omega²-2c²)r. These verify both force signs and
axis moments. An exact finite-amplitude translating entropy/rotation solution
has rho=f(z-Vt), Omega=Omega0/sqrt(rho), constant g=sqrt(mu0)c, Br=Bz=0 and
the same fixed p(r). All conservative variables are independently averaged;
it tests nonlinear transport of angular momentum and energy, not just a wave.

Small-amplitude torsional Alfven propagation is exact through first order;
nonlinear feedback is O(amplitude²) and is separated by amplitude controls.
Weak toroidal resistive decay uses J1(qr) exp(-Dq²t), with J0(qR)=0 giving the
conducting-wall trace. Its magnetic pressure/flow feedback is O(amplitude²).
These fixtures do not replace the original target or establish FRC formation.
