# Joint reaction and Coulomb momentum/energy source

`CollisionalAtomicMoments` combines elastic charged-species pairs and the resolved
atomic reactions. Every occupied charged heavy level
and the electron fluid retains independent momentum and lab translational energy.
Each unordered pair appears once. Neutral scattering and charge exchange are
absent, explicitly; no interaction is invented for zero charge or population.
The reference velocity is fixed as in atomic_frame.md.

The collision physics is the displaced isotropic Maxwellian Landau moments in
maxwellian_pair_collisions.md, based on
[Pfefferlé, Hirvijoki & Lingam (2017)](https://arxiv.org/abs/1701.08037).
The algebra below is a reformulation for conservative extensive energies, not a
different scattering law. For S=Ta/ma+Tb/mb (T in joules), u=va-vb, M=ma+mb,
mu=ma*mb/M, C=4π na nb (qa qb/(4π epsilon0))² lnΛ/mu:

    K = C [erf(x)-2x exp(-x²)/sqrt(π)] / |u|³, x=|u|/sqrt(2S)
    H = C sqrt(2/π) exp(-x²) / (M S^(3/2))
    R_a = -K u, R_b = -R_a
    dE_a/dt = R_a·Vcm + H (Tb-Ta), dE_b/dt = -dE_a/dt
    Vcm = (ma va+mb vb)/M.

Thus thermal sources are Qa=(mb/M)K|u|²+H(Tb-Ta) and
Qb=(ma/M)K|u|²-H(Tb-Ta). These equal the existing Gaussian moment formulas:
use B-|u|² A=sqrt(2/π) exp(-x²)/sqrt(S). Pair energy and momentum sums vanish
without an energy correction. The difference in bound internal levels does not
change during this elastic channel; inelastic reactions remain separate.

At each predicted midpoint, freeze K/(na nb), H/(na nb), velocities and
populations to construct a linear operator whose value at that state equals
these nonlinear sources. For normalized numbers f, scaled momentum p and energy
E in eV, substitute T_a=2[E_a-(p0/(2e))v_a·p_a]/(3f_a).
Momentum rates are -nu_a p_a+nu_b p_b, nu_a=K0 N f_b/ma.
Energy exchange rates use lambda_a=2 H0 N f_b/3 and opposite pair rows.
No division by tiny own population or kilogram mass-density is required.
All charged pairs advance together by exponential midpoint while their numbers
are constant. The full source uses a symmetric collision-half / reaction-full /
collision-half step. Reaction heat uses the original full-step midpoint. No
moment is reset between sources. Both subsolvers and the full split are second
order; this does not guarantee positivity at arbitrary dt.

The initial unsplit frozen-generator attempt passed the independent collision
limit but rejected the first reacting argon step: newly born rare charge states
had negative thermal energies. That source/failure is archived in outputs/
atomic_collisions. It is not repaired. Symmetric source splitting avoids freezing
collision temperature coordinates across population creation; its full-history
and time-halving evidence must be checked before adoption. `generator` still
returns the actual summed instantaneous source for diagnostics; `step` uses the
explicit split, rather than claiming its exponential is an exact coupled solve.

The analytic zero-drift limit and the existing small-x series evaluate drag.
One zero temperature is allowed when the combined variance is positive; both
cold with nonzero drift use the finite cold-beam moment limit. Both cold with
zero drift reject as singular. This is not a density/temperature floor or a claim
of weak-scattering validity for cold matter. Cold occupied species are reported.
Fixed positive lnΛ is supplied explicitly; no screening, magnetization, neutral
cross section, anisotropy, photon escape or physical wall model is added.

Linear independent mass/charge/momentum/energy invariants remain those of the
reaction model. Populations, current and temperatures change self-consistently
within this local source. `VectorizedCollisionalAtomicMoments` now connects it
to the precise moving planar reaction runner; the driven coil path still needs
changing-population integration. `add_pair_generator` exposes the unchanged
opposite-row assembly for density-free K0,H0, allowing the explicit elastic
neutral wrapper to use the same conservative coordinate representation. The
base class itself still selects only occupied charged pairs.
