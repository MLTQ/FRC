# Collisional magnetic decay with evolving resistance

Evolves the complete moving SpitzerWallPlasma for12us in a filled periodic
cylinder: R0.035m,L0.24m,Ar+,n1e21/m³,both gamma5/3,initial Te5/Ti1eV.
Velocity initially zero; Bz=A J0(qr),psi=A r J1(qr)/q,J1(qR)=0,zero bias.
The EOS entropy encoding includes its reference-density offset explicitly.
Default A=1e-4T is a small perturbation, with weak magnetization. Uniform
background electrons cool into ions while resistance increases in time.

The independent small-amplitude reference uses the separated scalar Coulomb
temperature clock, followed by a separate quadrature of electrical resistance:

    eta(T)=0.51 sqrt(me) e² lnLambda /
           [3 (2 pi)^(3/2) epsilon0² (e Te)^(3/2)]
    Bz,Jtheta = initial profiles * exp[-q² integral eta dt/(mu0(1+de² q²))]
    de²=mu/(mu0 n e²), mu=me mi/(me+mi).

No production collision or transport routine produces the reference. Its
temperature branch is restricted to1–5eV in this example. Nonlinear Lorentz
forces, motion, current stress and Joule heat remain active in the simulation;
their finite-amplitude corrections are absent from the linear reference.
The density is higher than the earlier constant-resistance fixture to keep
thermal relaxation relevant on the magnetic decay time and reduce mean free
paths relative to1/q. It is a mathematical component experiment, not a choice
of laboratory fill or an ionization-equilibrium assumption.

Weak-field scalar resistance, leading-mass-ratio electron Joule heat and
Maxwellian exchange are conditional closures. Conduction/thermal forces,
frictional heat flux, viscosity and neutral physics remain absent. Uniform
zeroth-order temperature and an axisymmetric azimuthal current make those
omissions vanish from the reference's first-order magnetic equation; nonlinear
heat redistribution is not a complete Braginskii prediction.

Independent order20 physical reconstruction measures B/J relative to both
initial and remaining signals, thermal histories, and physical energy. Energy
drift is scaled by the initial magnetic-plus-current energy, so a large thermal
background cannot hide it. Native energy and instantaneous power are separate.
Actual Hall/drift/log/rate diagnostics are saved at every endpoint.

Existing nr8/16 field grids, nr/2 fluid elements, degree5 and axial16 are used.
GuardedModalGauss keeps original1e-12 absolute/1e-10 relative bounds. Periodic z,
impermeable wall, fixed psi and zero Btheta boundaries are unchanged. No added
numerical damping, floor, field filter or adaptive retry. All endpoints/stages
and the source closure are retained. Own4GiB RSS/256MiB disk checks are discrete
after setup/steps. No unit tests, new verification framework or automatic replay.

Original8/80 and16/80 runs complete. Original8/160 rejects step17 on floating
wall-trace cancellation; the failure and its targeted observational replay
are retained. Use examples/lifted_spitzer_decay for the same experiment with
coordinate addition before magnetic reconstruction; its80/160 full runs
complete with unchanged bounds. See docs/SPITZER_DECAY_RESULTS.md.
