# Numerical thermal-energy transport in the moving pulse

## Higher-order central advance improves the actual hold (2026-10-16)

The localized central RK4 loss motivated an optional fixed twelve-stage DOP853
principal advance, implemented in `/private/tmp/frc-finite-drift` at `fbbaf8e9f`.
Main's running pulse source is untouched. The configuration option is
`central_integrator="dop853_fixed"` alongside the existing `energy_split` and
`central_energy_diffusion`. Omitted/`rk4` preserves the prior method. The complete
physical derivative is evaluated at every stage: momentum, coordinates, rotation,
both entropies, finite-mass current, prescribed coils, collisions and vacuum.
No energy correction, clipping, atmosphere or changed physical coefficient is
introduced. Read the private module's `plasma_bench/material_rk8.md`.

The tableau is supplied by installed SciPy 1.18.1's DOP853 class. Its principal
method is eighth order for smooth ODEs; see the
[official SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.DOP853.html).
We use fixed input steps, without embedded error control or dense output.
The unchanged outer diffusion splitting remains at most second order; this
is not an eighth-order claim for the whole solver or for nonsmooth terms.

Run 94170 exited 0 in 296.846 s. It reuses the identical actual step-1450
initial arrays, physical inputs and 43.75 ns constant-current hold. Both complete
50/100-step trajectories, exact source, arrays and inspected PNG/PDF are in
`outputs/material_hold_interval_dop853`, source
`df92f8cb55d196b0e89b3fad787044c4027be7585411ef99ee39c245bdfaab92`.

| Central method | Outer step | Energy change over hold |
|---|---:|---:|
| Prior RK4 | 0.875 ns | -1.1661176663e-4 J |
| Fixed DOP853 | 0.875 ns | -2.7961234994e-8 J |
| Fixed DOP853 | 0.4375 ns | -4.5601994403e-10 J |

The same-step energy loss falls by about 4170 times. Directly measured diffusion
energy increments are +8.94e-15/+1.62e-14 J, leaving the small remaining loss
in the central advance. Against the finer higher-order trajectory, the coarse
endpoint has relative B/J/velocity differences 1.42660e-5/6.30950e-5/7.89530e-6
and maximum material-position difference 7.69328e-9 m (6400 physical quadrature
points, actual finer volume). These differences are small in this particular
hold; they do not establish whole-pulse, mesh or experimental accuracy. Earlier
RK4 results and failed methods remain preserved, not overwritten.

Two complete configured pulses are now LIVE in the private checkout:
- 91606: weak 0.5 mT, 4x8, 800 steps, 4.375 ns;
  `outputs/material_tube_weak_rk8_energy_diffusion_r4_z8_t800`.
- 86765: strong 20 mT, 4x16, 4000 steps, 0.875 ns;
  `outputs/material_tube_independent_20mt_rk8_energy_diffusion_r4_z16_t4000`.
Their imported source is
`f5ddb11815a60a11f10dd650fb61180e305f22df7bcc44ee3a84ea1c87a428f2`.
Only the central integrator input changes from their corresponding prior runs.
The old-method finer full pulse 82391 remains LIVE in MAIN, source577ad923...bd0f.
Keep both source closures frozen. Collect real terminal results, then use the
existing configured-field comparison, which now names central_integrator too.
Do not confuse a central-method comparison with timestep convergence. All
physical heat/reaction/radiation/edge and formation acceptance gates stay open.

## Completed strong pulse: local loss is in central time integration (2026-10-16)

Strong pulse 41479 exited 0 after all 4000 steps through 3.5 us in 3524.619 s.
Its complete states and source 577ad923bec574d77ff1d2fb329047eba8cd7be539acdb26090603748343bd0f
are preserved in `outputs/material_tube_independent_20mt_energy_diffusion_r4_z16_t4000`.
Final external-work residual / E0 is -0.002153588 (maximum absolute 0.002154535).
Collision heat is 88.774993 mJ; numerical bulk ion heat is 33.281814 mJ.
Final mean Te/Ti are 9.745415/2.063095 eV, particle RMS radius 22.838098 mm,
middle/maximum boundary radii 55.988203/64.654447 mm. Completing the trajectory
does not establish temporal, spatial or physical accuracy.

Run 72041 exited 0 in 386.734 s. The narrow existing-state experiment
`examples/material_hold_interval.py` starts all three full continuations from
the actual step-1450 arrays and evolves the same constant-current interval,
1.26875 to 1.3125 us. Both prescribed source-current derivatives are zero, so
the independent energy limit is constant fluid + relative-current + induced-field
energy. All physical forces, collisions, numerical bulk viscosity, exact exterior
vacuum and numerical thermal mixing remain. There is no state or energy repair.
All accepted arrays, operator energy increments, source archive and inspected
PNG/PDF are in `outputs/material_hold_interval`, source
afe4723d925c814341b518a8cefc4aaaba6653f80bf33441daaa12ed9c19c92d.

| Solver step | Full steps | Total energy change | Diffusion contribution |
|---|---:|---:|---:|
| 0.875 ns | 50 | -0.116611767 mJ | +1.1713e-14 J |
| 0.4375 ns | 100 | -0.006434984 mJ | +2.5702e-14 J |
| 0.21875 ns | 200 | -0.000243639 mJ | +4.0024e-14 J |

The energy loss occurs during the central RK4 steps, and decreases by factors
18.12 and 26.41 under halving. This local result identifies finite-step central
integration as the source of the observed hold loss; it does not prove a
particular asymptotic order or exclude other full-pulse errors. The coarse
endpoint reproduces archived energy within 6.0e-14 J and positions within
2.2e-14 m; tiny timestep arithmetic differences prevent bitwise identity.
Independent centered energy differences along the actual central derivative at
three saved states decrease approximately quadratically with perturbation size
before roundoff, reaching <=8.3e-4 W, far below the coarse average loss ~2665 W.

Actual endpoint fields use 6400 physical quadrature points and the finest
continuation's volume. Relative B/J/velocity differences against 0.21875 ns are
6.0401%/33.1348%/4.8805% for 0.875 ns, versus
0.190964%/1.016104%/0.108900% for 0.4375 ns. Maximum position differences are
89.1465 and 1.24069 micrometres. These local improvements do not establish
full-pulse convergence, particularly through later evolving fields.

Existing numerical-method comparison 99873 also exited 0, using 81 common
times and 6400 physical points. Output
`outputs/material_tube_energy_diffusion_vs_upwind_r4_z16_t4000` preserves JSON
and inspected PNG/PDF. Against the previous upwind RK4 run, final relative
B/J/velocity differences are 55.3026%/92.4714%/12.6432%, with maximum material
position difference 7.64067 mm. Both the numerical spatial flux and integrator
change, so these are method sensitivities, not timestep convergence or a ranking
against truth. The new method cannot yet be claimed more accurate.

Full-pulse 82391 remains LIVE in MAIN at 8000 steps requested (0.4375 ns),
unchanged moving source 577ad923...bd0f. Keep its imported source frozen and
collect its actual terminal result before the existing configured-field
comparison. The private checkout is now idle. If improvement is still needed,
use the localized central time-integration error to motivate a concrete solver
change; do not substitute repeated blanket refinements or more report tooling.
Physical heat, reactions, radiation, edge conditions and formation gates remain.

## Strong-pulse conservation observation and matching time control (2026-10-16)

The weak-case agreement below does not extend to demonstrated strong-pulse
accuracy. Live handle 41479 reached 3500/4000 steps (3.0625 us). Its maximum
observed absolute external-work residual / E0 is 0.0021545345 at 2.275 us,
with E0=0.0895071115121 J. Most observed loss occurs near the start of the
constant-current hold. Between saved times 1.26875 and 1.3125 us, the bias/array
vacuum-center fields stay [0.0001,-0.02] T, integrated external work stays
0.253495534314572 J, but total energy falls from 0.342987094621427 to
0.342870482854735 J: a 0.116611767 mJ loss with zero source-work increment.
This interval cannot be explained solely by quadrature of a changing source.

The added mixing halves' cumulative represented thermal-energy residuals remain
around 1e-14 J per species, far smaller than the full-system loss. The dual
entropy sampling matrix is one-hot, so its averaged h energy identity does apply
to represented entropy. These observations do not by themselves identify the
entire numerical cause or prove coupled physical accuracy.

Matching full-pulse handle 82391 is now running in main, changing only 4000 to
8000 steps (0.875 to 0.4375 ns). Its new config is
`configs/material_tube_independent_20mt_energy_diffusion_t8000.json`; output is
`outputs/material_tube_independent_20mt_energy_diffusion_r4_z16_t8000`.
Both runs use unchanged source
`577ad923bec574d77ff1d2fb329047eba8cd7be539acdb26090603748343bd0f`.
Keep the main and private moving source closures frozen while they run. Collect
actual terminal outcomes and use the existing field comparison. This control is
motivated by measured strong-pulse energy loss; no convergence claim follows
from its launch or from the small weak-pulse differences.

## Completed weak pulse and terminal pressure-switch comparison (2026-10-16)

The implementation is now integrated into main as `f838aeab0` (the unchanged
local implementation commit was `e0f7c5376`). The private checkout remains frozen
while strong-pulse handle 41479 runs. Earlier running-status paragraphs below
are historical; only 41479 remains live at this checkpoint.

Both complete 0.5 mT, 4x8 weak-array controls finished through 3.5 us: handle
9903 exited 0 after 800 steps in 470.061 s; 17128 exited 0 after 1600 steps in
878.699 s. Their frozen source is
`577ad923bec574d77ff1d2fb329047eba8cd7be539acdb26090603748343bd0f`.
Full histories are in `outputs/material_tube_weak_energy_diffusion_r4_z8_t800`
and `outputs/material_tube_weak_energy_diffusion_r4_z8_t1600`.

Existing configured-field comparison 38599 exited 0. Its inspected PNG, PDF and
JSON are in `outputs/material_tube_weak_energy_diffusion_timestep_comparison`.
Only dt changes: 4.375 to 2.1875 ns. The comparison uses 81 common actual times
and 3200 material quadrature points, with no temporal interpolation or repairs.

| Measured quantity | 800 steps | 1600 steps |
|---|---:|---:|
| Maximum observed external-work residual / E0 |3.07170e-8|1.01524e-8|
| Maximum observed response-work residual / E0 |7.10122e-8|2.27172e-8|
| Final collision heat |0.154041857 mJ|0.154041782 mJ|
| Final particle RMS radius |27.492161008 mm|27.492161010 mm|
| Final maximum/middle radius |54.424296608 mm|54.424296561 mm|
| Final mean Te / Ti |3.493268384 / 0.375930388 eV|3.493268382 / 0.375930388 eV|
| Numerical bulk-viscous heat |0 J|0 J|

Final relative vector differences are total B 5.63726e-6, induced B 1.26301e-5,
J 1.11669e-5 and velocity 4.43923e-9. Maximum material-position mismatch is
1.88272e-10 m. Across the sampled pulse, maximum total-B difference is
7.98856e-6, induced-B difference 0.00104722 and J difference 0.00119420;
all three maxima occur at 0.2625 us. Small early reference signals affect the
component-relative ratios. These are small measured time-halving differences
for this weak case, not a spatial-convergence or strong-pulse result.

Fine-run cumulative numerical energy redistribution is
[9.0517794e-6, 1.3166376e-11] J for electrons/ions. This measures redistribution,
not generated heat. Corresponding entropy production / kB is
[2.7122920e13, 4.2180027e8], with represented thermal-energy residuals
[-6.96610e-15, -1.22513e-15] J. Physical nonlocal heat remains absent.

The old original-resistance pressure-switch run 95749 exited 1 after 3773
accepted steps, at 1.6506875 us. Attempt 3774 failed with
`Upwind midpoint negative-pressure branch did not converge`. Its accepted prefix
and failure are preserved without repair; see MATERIAL_COMPRESSION.md. The
strong 20 mT energy-diffusion run 41479 is still running. Collect that actual
terminal result before interpreting this method's strong-pulse behavior.

The third timestep refinement did not settle detailed fields: paired final B/J
differences increased from10.19%/20.91% to12.44%/25.36%. The existing upstream
entropy flux is continuous through zero current, but its energy-paired pressure
force can jump there. Integrating that switch implicitly is already under study
in the original-resistance run95749. It is not yet a proven explanation for all
observed spatial and temporal sensitivity.

## Solver change and its scope

An alternative numerical spatial flux is implemented on local branch
`codex/finite-drift-collisions`, commit `e0f7c5376`, in
`/private/tmp/frc-finite-drift`. It retains the central species entropy/pressure
work adjoint and every physical current, fluid, collision and prescribed-source
term. The added diffusion transfers internal energy between material dual
volumes without exchanging energy with the magnetic/current state. There is no
sign-dependent numerical pressure force. This changes the numerical scheme;
it is not a new physical heat conductivity or a way to repair total energy.

Let m_i be fixed actual dual-volume particles, h_i internal energy per particle,
q the oriented relative particle flux, c_s the electron/ion drift fractions and
I the edge incidence. With fixed geometry, h_i=a_i exp(s_i). The added equation is

    M dh/dt = -I.T diag(|c_s q|/2) I h.

It conserves sum(m_i*h_i). For each species, the physical entropy production
divided by kB is the sum over edges of

    |c_s q|/2 * (h_head-h_tail)^2 /
               [(gamma-1)*h_head*h_tail] >= 0.

The Metzler generator's exponential preserves positive internal energy. The
code uses its sparse exponential action, then the EOS relation log(h_new/h_old)
to update entropy. It rejects nonpositive/nonfinite results without clipping.
The reported energy residual uses the actual represented entropy, never a target
normalization. Local entropy minima need not be preserved under heat transfer.

On a uniform straight mesh the diffusion coefficient is |u_s|*dx/2, vanishing
with grid spacing for smooth fields. This is declared first-order numerical
diffusion along relative-current edges, not physical anisotropic/nonlocal heat
conduction. Both numerical spatial error and full coupled time error must be
measured; positivity and energy conservation alone do not prove accuracy.

Diffusion half-steps surround full central RK4; smooth-time order is at most two.
Geometry/current are fixed only during each diffusion half. Actual collision
drag and both species heats remain in the central stages. Exact exterior vacuum,
the moving material inventory, and the existing wall-contact rejection remain.
No escaping-particle, sheath, neutral, ionization or radiation physics is added.

The config requires both `integrator="energy_split"` and
`entropy_transport="central_energy_diffusion"`; older inputs retain their
original solver. The existing field comparison now names both numerical inputs.
Both independent source-work budgets remain defined for the central pressure
law. Numerical thermal redistribution, entropy production and residual are
recorded separately from actual collision heat and numerical bulk-viscous heat.
Full derivation/contracts are in the new module companion in the private checkout.

## Completed independent mathematical evolution

Run93705 exited0 in0.891s, source
`84047d488b6f16f9c86bd37aeea54653a3d3dcedddb5400319218568dec8f3a8`.
`outputs/material_energy_transport` preserves all51-time histories, analytic
references, source archive and inspected scientific PNG. This is the added
diffusion component alone, not an assertion that the whole plasma is verified.

A unit periodic domain starts with exact cell averages of h=1+0.8*cos(2*pi*x).
At fixed D=0.1 the independent continuum solution decays as exp(-4*pi^2*D*t).
Final L2 errors at t=0.5 for16/32/64/128 cells are1.99547e-3,4.98472e-4,
1.24592e-4,3.11464e-5, a second-order spatial decrease. Maximum relative
thermal-energy drift is2.73115e-14; sampled total-entropy increments are positive.

With fixed |u|=1 and D_num=dx/2, the change from the initial profile decreases
0.25635,0.14948,0.08075,0.04197 as numerical diffusion vanishes. An independent
unequal-particle two-reservoir exponential solution agrees to5.32907e-15
absolute energy per particle; its relative total-energy change is1.55431e-15.
No fitted coefficient, state repair or unit-test suite was used.

## Coupled trajectories still running

All new pulse runs use frozen source
`577ad923bec574d77ff1d2fb329047eba8cd7be539acdb26090603748343bd0f`
in the private checkout. Do not change its imported source until they terminate.

-41479: full20mT paired-collision array,4x16,4000 steps/3.5us,0.875ns;
 outputs/material_tube_independent_20mt_energy_diffusion_r4_z16_t4000.
-9903 and17128: complete0.5mT weak-array controls,4x8,800/1600 steps over3.5us;
 outputs/material_tube_weak_energy_diffusion_r4_z8_t800 and t1600.
 Physical inputs match between these controls; only the timestep changes.
-95749: the older original-resistance pressure-switch split+bulk-pressure run
 remains live in MAIN with unchanged source31721f1f...d17afc.

An initial3x6 weak-control configuration was rejected before evolution because
the underlying mesh requires at least four cells per direction. Its terminal
failure, exact input and source archive are preserved in
outputs/material_tube_weak_energy_diffusion_r3_z6_t800. It is not a physical
failure or a completed control; the new4x8 runs have separate output directories.

Collect the exact live handles. After both weak runs terminate successfully,
use existing `examples.material_coil_plot --configured-fields` in the private
checkout to measure the complete800/1600 comparison. Collect the full strong
pulse before deciding whether this numerical change improves field accuracy.
Keepfrc-vacuum.3.4/.3.3, transport and the entire original FRC goal OPEN.
