# Shared smooth finite-mass transport and electric-kernel representation

This selectable reference addresses two errors exposed by the conservative
moving baseline: inaccurate near-solenoidal mass/enthalpy transport and nonlinear
electron entropy advection. It retains the parent finite-mass equations, physical
masses, full energy/current terms, positive-density periodic cylinder and explicit
RK4 of FINITE_MASS_NUMERICS. It does not replace the original source or enable
the formation workflow. Its magnetic cancellation error remains unresolved.

## Conservative reconstruction and face work

Use x=r²/2. Actual annular cell averages are fitted to r^s P3(x), with s=0 for
even scalar quantities and axial momentum, s=1 for radial momentum and Br, and
s=2 for angular momentum. The mass-weighted moment matrices are integrated
exactly. Shared face interpolation uses four physical interior cells; periodic
axial interpolation uses four cell averages. No exterior matter is invented.

Nonlinear face products use physical Gauss quadrature. Transverse cell
polynomials average two mirrored cubic stencils, each of which includes the
target cell and integrates back to its original average. This distinction is
essential: switching an arbitrary-point face stencil halfway through a cell
does not preserve its transverse integral. Independent arbitrary-data and
polynomial tests check mean preservation through the axis and wall.

At each face quadrature point, one mass flux F_rho=m_normal supplies

    F_mr = vr F_rho - Br Bnormal/mu0 + delta_normal,r (p+B²/(2mu0)),
    F_mz = vz F_rho - Bz Bnormal/mu0 + delta_normal,z (p+B²/(2mu0)),
    F_L  = (L/rho) F_rho - r Btheta Bnormal/mu0,
    F_W,bulk = [(Kb+Kj+Hi+He)/rho] F_rho.

The same pressure, rotation and magnetic reconstruction supplies the cylindrical
radial stress source. Full computed Poynting power and mixed current energy,
enthalpy and drift-stress fluxes are added by the existing conservative contract.
The actual staggered normal B face averages are retained. Corner flux still
defines a divergence-free field; high-order mechanical traces do not change
the older electromagnetic degrees of freedom or energy Hessian.

Electron and ion pressures are reconstructed separately from their existing
dual/primal representations. Current energy is interpolated as a face quantity;
tiny signed interpolation overshoots are not a stored physical state or a claimed
positive current-energy quadrature. Stored current energy remains its positive
physical J²/n form. Density and both reconstructed pressures must be positive,
or the stage rejects. Neither cell energy nor thermal state is repaired.

This is a **central smooth reference** with no Rusanov dissipation, viscosity,
resistivity, trace limiter or hidden fallback. It has no shock admissibility,
entropy inequality or general nonlinear stability theorem. The unchanged
limited solver remains available and its measurements retain their own source.

## Shared electron entropy transport

The new primal mass flux still provides the only continuity equation. Its exact
overlap supplies dual n_t and particle flux Gamma. Reconstruct conservative
dual n and U separately, then use

    U_t = -div[Gamma (U_face/n_face)] + current transport + heat source.

Reconstructing U/n directly as a cell average introduces a nonlinear averaging
error. On nonuniform radial dual x-volumes, cubic face interpolation matches
their true volume moments. Axial interpolation is periodic. Constant U/n follows
any shared particle flux exactly. For uniform axial translation, the primal
flux map commutes with interpolation and Gamma/n_face is the common velocity,
so both n and U obey the same high-order linear advection operator.

The exact contact has n(z-Vt), pe(z-Vt), pi=p0-pe and uniform V/Bz. Its full
Ez=-pe_z/(ne)=pi_z/(ne) balances the separate species forces; it is curl-free,
possibly with a periodic harmonic voltage. No forcing or heat is needed.
Initial nodal pe is analytic, while cell pi complements the actual Q1 electron
average to keep total pressure constant. The resulting initial pi/energy
projection errors are measured separately. Evolution never resets this partition.

An explicit limited-entropy option is retained only as a controlled comparison.
Neither smooth nor limited advection has a proven physical entropy inequality.
Nonpositive reconstructed n/U rejects, without floors or changing the mass flux.

## Exact electric-kernel split and honest residuals

The original contact rate fails the old relative induction residual check even
though its full-equation backward error is near roundoff. The electric field
has norm about 61.4 in coefficient units, while its non-kernel remainder is only
2.20e-10. Recombining those components in double precision cannot retain a
1e-10 relative accuracy for the much smaller curl. The original failed source,
input, forcing and exception are preserved.

The new solver algebraically separates Er/r independent of z and Ez independent
of r. Both are exact kernel fields of the stated Faraday operator. There is no
amplitude threshold, removal of a physical signal or changed electron mass.
Differences are taken before division by metric factors:

    psi_t = -r Etheta,
    g_t = (Er_coeff[j]-Er_coeff[j+1])/dz
          +(Ez[i+1]-Ez[i])/Delta x.

The extracted fields vanish bit for bit under this form. Solve the non-kernel
problem with every original strict relative Ohm/Ampere/Faraday/operator check
at the unchanged tolerance. Restore the complete E for Poynting transport and
additionally check the original full equations by backward error, using their
un-cancelled term magnitudes. Report the full-array Faraday error relative to
the tiny rate separately: it is about 1.67% in the preserved contact example,
while original backward errors are below 9e-17. It is not described as a tiny
rate-relative residual. Weak physical signals beside a large kernel are tested
against an independently assembled original mixed system, with their input
roundoff amplification explicitly bounded.

This separable extraction is not a general Helmholtz decomposition. Mass inversion
and subtraction still produce roundoff-sized input noise, and more general nearly
curl-free fields may reject. Failures remain atomic. The old moving-current and
inertia solver sources remain unchanged; only the smooth reference selects this
representation. All reported residual and conditioning components stay distinct.

## Why the slow magnetic mode still fails

Improved pressure transport does not establish an accurate magnetic rate. In the
slow linear mode, -v×B0 and delta J×B0/(ne) are each about 19200 times their net
electric field. The existing mass-weighted nodal velocity map and weak current
representation do not reproduce that cancellation accurately. At nr64, the
magnetic-rate error remains about 16.84 times its physical signal.

Independent linear controls retain the original n, physical mass, B0 and inertia
but substitute exact velocity or current only inside the forcing. Exact velocity
alone reduces nr64 error to 3.21; both exact give .01991. These are diagnostic
oracles, **not** acceptable production patches or changed physical initial states.
The rate-only linear control uses perturbation-only magnetic input to exclude
bias-current roundoff; production trajectories retain their full original fields.

The next correction must make velocity/current/induction representation consistent
and more accurate, with complete physical B_t/J_t/E norms and local energy tests.
A higher-order compatible magnetic/electric representation is a candidate to
derive and verify. Do not initialize velocity from numerical current to hide
the error, insert exact fields into evolution, or judge success from pressure
and global conservation alone.

Physical collisions, heating partition, conduction/ionization, the original
heated free edge and exact vacuum, real ends/spatial coils, shock and kinetic
validity, converged formation, three-dimensional stability/confinement and the
experimental interface remain separate unfinished requirements.
