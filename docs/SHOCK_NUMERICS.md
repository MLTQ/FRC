# Magnetic shocks, numerical heating and implicit-solve accuracy

The fixed-chamber coupled solver's centered magnetic update has no inherent
shock diffusion. Its exact discrete work pairing does **not** establish a
correct entropy solution or a conservative local magnetic momentum flux. An
independent strong perpendicular shock exposes magnetic oscillations and even
a false field reversal while the global energy balance remains accurate.
Those fields cannot be interpreted as FRC formation or physical reconnection.

## Independent cylindrical standing shock

An inward radial flow with axial B has a perpendicular shock at r_s=0.5 m.
Upstream at the shock, rho=1 kg/m³, p=0.1 Pa, v_r=-6 m/s, B_z=sqrt(mu0) T,
gamma=5/3. Physical eta=0. The downstream state is the unique compressive,
entropy-increasing root of mass, momentum, induction and total-energy jumps:

```
[rho v] = 0
[rho v² + p + B²/(2 mu0)] = 0
[v B] = 0
[v (gamma p/(gamma-1) + rho v²/2 + B²/mu0)] = 0
```

The root is independently constructed in `shock_reference.py`; it is not
extracted from either numerical solver. The perpendicular jump conditions and
admissible compression are described in
[Fitzpatrick's perpendicular-shock derivation](https://farside.ph.utexas.edu/teaching/plasma/Plasma/node104.html).
Our cylindrical continuation follows directly from the steady ideal-MHD equations:

```
r rho v = m; B/rho = k; p/rho^gamma = K (separate K on each side)
H = gamma K rho^(gamma-1)/(gamma-1) + m²/(2 r² rho²) + k² rho/mu0
d[r(rho v²+p+B²/(2mu0))]/dr = p+B²/(2mu0)
```

The scalar Bernoulli equation selects the sub-fast downstream or super-fast
upstream branch. It rejects radii past choking. Finite-difference differentiation
of the independent branches checks the cylindrical momentum equation; a separate
zero-field test reproduces the analytic hydrodynamic compression ratio.

The exact branches occupy 0.38<=r<=0.85 m in an R=Z=1 m chamber. An inner
extension gives regular radial velocity at the axis. An outer cosine taper makes
velocity zero at the reflecting wall. These extensions are not exact steady
solutions. Comparisons use only 0.44<=r<=0.66 m, |z|<=0.1 m through t=0.02 s.
Independent controls move the inner extension, double chamber length and remove
the remote wall taper in a shorter t=0.01 s run. The last control is shorter
because unmodified inward supersonic flow at a closed outer wall produces a
separate difficult rarefaction. No free-vacuum robustness claim follows.

Initialization uses 12-point Gauss integration of rho, rho v, gas energy and
r B. The jump aligns with a radial face. Density and B references are cell
volume averages. Pressure comparisons recover p from the averaged conserved
gas variables, including their O(h²) subcell kinetic-energy variance. Error
weights use actual cylindrical overlap volume of the declared window.
The bounded conservative reference has a conducting radial wall and periodic
ends; the coupled solver has reflecting gas boundaries and an infinite magnetic
exterior. Neither is a global steady-flow boundary-value solution.

## Explicit numerical magnetic diffusion

The optional dimensionless coefficient `shock_diffusion=C` adds

```
D_num = C h min(h max(-div(v), 0), |v| + sqrt(gamma p/rho+B²/(mu0 rho)))
eta_num = mu0 D_num
```

Here h=min(dr,dz). A centered cylindrical velocity divergence uses volume
centroids; linear radial flow gives its analytic regular-axis limit. Nonnegative
adjacent-cell averaging maps eta to plasma current nodes. Vacuum and exterior
modes have precisely zero resistivity and no fluid. Diffusivity is O(h²) for
smooth compression, O(h) at a compressive shock, and zero for rest/expansion.
The sensor uses midpoint velocity/field iterates and entering gas pressure in
its speed cap. Numerical diffusion remains **off by default**, is not a
conductivity closure, and does not model physical reconnection or transport.

It enters the same positive nodal resistance as physical eta. Physical and
numerical heat are accumulated separately:

```
Q_phys/num = dt (2 pi/mu0) sum[R_phys/num (K u_mid)²] >= 0
Delta W_induced + Delta K_gas + Q_phys + Q_num = prescribed-source work
```

Both heat terms go into gas internal energy through nonnegative volume shares.
This is magnetic energy loss deposited as heat, not a posteriori pressure or
energy repair. Global Q_num includes remote extensions and cannot be labeled
shock-only dissipation. Even as a coefficient vanishes with h, integrated
dissipation at a shock need not vanish; an entropy-producing weak limit can
have finite dissipation. Profile/jump convergence is still required.

Coefficient, time and mesh variations are reported. The sensor does not guarantee
monotonicity or handle every MHD discontinuity. Finite-grid B/rho errors and field
overshoot remain measured limitations. The centered gas/magnetic split has no
new proof of general local momentum conservation or correct shock limits.

## Faster solution of the same midpoint equations

The retained `midpoint_method="direct"` factors the velocity-eliminated sparse
linear system at each nonlinear iteration. The default `"picard"` evaluates
the same implicit equation as a fixed point, then completes exterior nodes and
modes with a cached vacuum factorization. At each candidate it recomputes
advection, midpoint velocity, numerical resistance and the nonlinear residual.
Only converged iterates meeting change/residual tolerances may be accepted.
Default relative tolerance is 1e-13; failure leaves the accepted state untouched.
The time cap includes physical plus estimated numerical resistivity as well as
wave speeds. It is an accuracy/convergence restriction despite implicit time
integration. There is no automatic acceptance of a failed Picard iteration.

Direct/Picard comparisons cover both magnetic boundaries and independent
physical/numerical resistance choices. Their equations and converged work
identity are unchanged. Recorded runtimes are CPU measurements on the stated
machine, not a claim of GPU performance or general real-time operation.

See [the measured shock report](SHOCK_VERIFICATION.md). Original Stage 2
formation remains gated; moving plasma interfaces, transport/ionization,
appropriate end conditions and regime-specific Hall/kinetic validity remain
required in addition to stronger numerical verification.
