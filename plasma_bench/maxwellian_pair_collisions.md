# Finite-drift Coulomb momentum and thermal moments

`pair_moments` returns paired interspecies force and heat for positive-density,
isotropic displaced Maxwellians, nonrelativistic velocities and weak small-angle
scattering. It takes an explicit positive Coulomb logarithm; no cutoff model,
logarithm floor, neutral collision or vacuum atmosphere is invented. Arrays
broadcast over samples; the last drift dimension contains three components.

The kinetic basis is the Landau operator and Gaussian convolution in
[Pfefferlé, Hirvijoki & Lingam (2017), equations 11–24 and 40–48](https://arxiv.org/abs/1701.08037).
This implementation specializes to isotropic Maxwellians, not the paper's full
Hermite or ten-moment system. It is not a Spitzer/Braginskii conductivity tensor.

For temperatures Ta,Tb in joules, define M=ma+mb, mu=ma*mb/M,
S=Ta/ma+Tb/mb, U=Ua-Ub, x=|U|/sqrt(2S) and
C=4*pi*na*nb*(qa*qb/(4*pi*epsilon0))²*lnLambda/mu. Independent binary-scattering
reasoning gives relative velocity g and pair centre velocity G:

    R_a = -C <g/|g|³>, R_b = -R_a
    Q_a = -C <(G-Ua).g/|g|³>
    Q_b =  C <(G-Ub).g/|g|³>

The relative velocity is Gaussian with mean U and variance S. Conditional on g,
<G|g>=(ma*Ua+mb*Ub)/M + (Ta-Tb)/(M*S)*(g-U). Its integrals are

    A = <g/|g|³> = [erf(x)-2x exp(-x²)/sqrt(pi)] U/|U|³
    B = <1/|g|> = erf(x)/|U|
    P = C U.A = -R_a.U >= 0
    Q_a = [Ta/(ma*S)] P + C (Tb-Ta)/(M*S) B
    Q_b = [Tb/(mb*S)] P - C (Tb-Ta)/(M*S) B

Consequently total momentum is conserved and Qa+Qb=P restores lost relative
kinetic energy. Each heat can be negative during temperature equilibration;
neither is clipped. Total entropy source Qa/Ta+Qb/Tb is nonnegative. At zero
drift, A/U=sqrt(2/pi)/(3*S^(3/2)), B=sqrt(2/pi)/sqrt(S), recovering the existing
thermal equilibration coefficient. At large drift force decays as |U|^-2.

The code uses a small-x Taylor evaluation below x=0.05 for roundoff accuracy,
otherwise the regularized incomplete gamma function for the drag numerator.
This changes numerical evaluation, not the collision law. No state integration,
energy projection or stabilization is performed here.

`examples.maxwellian_pair_relaxation` exercises complete homogeneous collision
relaxation against direct velocity-space quadrature and conservation-derived
equilibrium. Existing moving-coil sources remain unchanged while live runs finish.
Later coupling must use one paired force/heat closure in every driven/undriven
Ohm solve and species source, replacing both prescribed drag and thermal exchange
without double counting. In the current solver w=Ui-Ue; use Ue-Ui=-w here.

Maxwellians need not remain Maxwellian under finite drift; stress anisotropy,
heat flux, kinetic instabilities, magnetized collision corrections and a justified
screening prescription remain open. Agreement of moments is not validity of a
five-moment closure for the formation pulse.
