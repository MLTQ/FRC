# Physical current transport and drift stress

This extension keeps the regular fluid and compatible magnetic/current fields
of REGULAR_COUPLING_NUMERICS. It adds the instantaneous inertial transport
terms, not a complete species-pressure or time-evolution model. The verification
domain is a filled cylinder with periodic z and impermeable radial wall. There
is no vacuum interface here; formation still requires exactly zero exterior
plasma. Physical electron/ion masses are unchanged.

Let w=J/(e n), reduced mass mu=mi me/(mi+me), delta=(mi-me)/(mi+me), and
T=mu n w w. The original terms are

    bulk force: f = -div(T)
    Ohm forcing: mu/e [(v·grad)w + (w·grad)v - delta (w·grad)w].

All gradients include the cylindrical connection: for axisymmetric vector a,
(grad a)_{r,theta}=-a_theta/r and (grad a)_{theta,theta}=a_r/r.
The axis has dependent regular transverse traces, not radius floors or extra
zero-energy velocity unknowns. The actual n_t is the original fluid continuity
rate, not its analytic target.

## Weak advection and consistency defects

For a current test field t, the volume form is

    A(t,J;v,n_t) = mu/(2e²) integral
       [(t·Dv J - J·Dv t)/n + (t·J) n_t/n²] dV.

At each field-cell interface, oriented from left to right, add

    mu/(2e²) integral a_hat [t_L·J_R - t_R·J_L] dS,
    a_hat = 0.5 [(v_L·normal)/n_L + (v_R·normal)/n_R].

The face form is also skew. Thus A(J,J) is exactly the required density-energy
exchange at the integration quadrature. Away from faces its strong field is

    mu/e (v·grad)w + mu J [n_t+div(nv)]/(2e² n²).

This last term is a measured local consistency defect. Skewness is not proof
of original-equation accuracy. Normal-current jumps and div(J) are separately
measured. Self drift is -delta A(t,J;w,0); its strong discrepancy is proportional
to div(J), and its current work is zero by the skew form. Neither defect is
silently discarded.

## Drift stress and stretching

The volume force is the physical velocity-gradient adjoint of integral
T:grad(v) dV using the actual declared velocity map. At radial interfaces use
shared traction 0.5*(T_L+T_R)·normal. Its bulk work subtracts
integral (v_L-v_R)·traction. Stretching includes the identical paired face
work, linear in the current test field with the actual J on the traction side.
Axial fluid velocity is continuous through its periodic Fourier reconstruction;
there is no corresponding velocity jump there. The symmetric stress gives
zero net axial force and angular torque in the closed periodic fixture.

Numerical radial transport flux is zero at the wall and axis. The analytic
fixture has zero normal current at the wall; reconstructed traces are measured
rather than assumed exact. The axis surface has zero cylindrical area. Radial
field cells subdivide fluid elements so quadrature never crosses an unaccounted
fluid discontinuity. No upwinding, viscosity, resistivity, density atmosphere,
filter, temperature reset or fitted global energy correction is added.

## Verification and physical limits

The study tests each original term independently against analytic fields with
complex derivatives. It freezes the calculation source, all coefficients and
metrics. Audit replays every array/metric, then uses independent root-product
polynomials and FFT reconstruction at orders 16 and 20. Force is also measured
between nodes using the exact derivative of the declared physical velocity.
Global, first/last-element and maximum errors are reported. Omitting faces is
an explicit negative accuracy control despite retained algebraic skewness.

Grid fluid energy and the exact continuous variable-density energy are not
identical. The audit therefore integrates the full reconstructed physical
energy derivative separately; it does not relabel grid conservation as exact
continuous conservation. Small instantaneous error for this fixture still
does not establish coupled trajectory stability or species-energy transport.
Species pressure/current-entropy, enthalpy/Poynting and unforced original-parent
wave/contact evolution remain adoption gates under frc-7mm.4 and frc-7mm.5.
