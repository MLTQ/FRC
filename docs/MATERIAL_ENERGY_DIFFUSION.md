# Numerical thermal-energy transport in the moving pulse

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
