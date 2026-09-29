# Infinite-cylinder analytic skin response

An independent continuum solution for a stationary cylinder of radius R with
uniform magnetic diffusivity D=eta/mu0. A prescribed exterior uniform axial
field is Re(Bs exp(i omega t)). The exterior has zero current. There is no
surface current: Bz is continuous at R; volume J can jump there. The solution
is axially invariant and initialized in its exact periodic state.

With k=(1+i)/delta and delta=sqrt(2D/omega),
`Bhat=Bs I0(kr)/I0(kR)`, `psihat=Bs r I1(kr)/(k I0(kR))`,
`Jhat=-Bs k I1(kr)/(mu0 I0(kR))`. Outside, psi extends as its boundary value
plus Bs*(r²-R²)/2, Bz=Bs, and J=0. Fields use exp(+i omega t).
Scaled Bessel functions avoid overflow in the interior ratios.

`mean_ohmic_power_per_length` integrates boundary Poynting work analytically
in time, independently of a numerical volume integral of eta J²/2. This is
stationary conductor verification, not a compressing plasma or a finite coil
experiment. The general magnetic-diffusion and skin-effect equations are in
MIT's *Electromagnetic Fields and Energy*, sections 10.5 and 10.7; the cylindrical
Bessel solution here follows by separation in the cylindrical diffusion PDE.
