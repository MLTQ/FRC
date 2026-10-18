# frc_scaling.py

## Intent
First-cut design relations for a field-reversed theta pinch (FRTP): pick a gas,
fill pressure and external field, get the equilibrium FRC, the formation and
pre-ionization requirements, and the equivalent single-turn coil/bank numbers.
It is explicit scaling, not a simulation. Its job is to choose the regime that
the benchmarked solver then refines.

## Equations (SI, temperatures in eV)
- Fill: diatomic gas, `n_fill = 2 p/(k T_gas)` atoms.
- Geometry: `r_s = (r_s/r_w) r_w`, `x_s = r_s/r_c`, `l_s = (l_s/l_c) l_c`,
  `V_s = shape * pi r_s^2 l_s`.
- Inventory: `<n> = f_N n_fill pi r_w^2 l_c / V_s`; `n_peak = <n>/<beta>`.
- Average beta (elongated FRC, flux-conserving coil): `<beta> = 1 - x_s^2/2`.
- Pressure balance: `<n_i> (T_i + Z T_e) e = <beta> B_e^2/(2 mu_0)`.
- Rigid rotor `B_z = B_e tanh(K(2r^2/r_s^2 - 1))`: trapped flux
  `phi_p = pi r_s^2 B_e ln cosh K/(2K)`.
- Bias needed: `|B_b| = phi_p/(f_phi pi r_w^2)`; the sign opposes the main field.
- Classical flux e-folding with a fixed profile:
  `tau_phi = mu_0 r_s^2 ln cosh K/(8 K^2 eta)`, from the Ohmic E_theta at the null.
  `eta` is NRL perpendicular Spitzer, `1.03e-4 Z lnL T_e^-1.5`, with the NRL
  Coulomb logarithm at `n_peak`. This is an upper-bound estimate; FRC
  experiments commonly lose flux faster than classical.
- Kinetic size: `rho_i = sqrt(2 m_i T_i e)/(Z e B_e)`, `s_bar = r_s/(4 rho_i)`
  (null-to-separatrix integral with rho_i at B_e). Conventions in the
  literature differ by O(1); `r_s/d_i` with `d_i = c/omega_pi` is also reported.
- Implosion: planar snowplow transit `t_sp = sqrt(2) r_w/v_A(B_e, rho_fill)`.
  Main-bank quarter period `t_q = rise_per_implosion * t_sp` (default 3).
- Coil: thin single-turn sheet, centre field `B = mu_0 I/sqrt(l^2 + 4a^2)`,
  Wheeler inductance `mu_0 pi a^2/(l + 0.9a)`.
- Bank (lossless LC into coil plus stray L): `E = L I^2/2`,
  `C = (2 t_q/pi)^2/L`, `V = I L pi/(2 t_q)`. Energy is set by the field and
  inductance; a faster rise costs voltage, not energy.
- Pre-ionization: the minimum T_e at which the column skin time
  `mu_0 r_w^2/eta` is `pi_skin_margin` (default 10) quarter periods.

## Empirical inputs (`Assumptions`)
Particle retention `f_N`, flux retention `f_phi`, `r_s/r_w`, `l_s/l_c`, the
shape factor, K and T_i/T_e are round mid-range values. They are not fitted to
any device. The benchmark (beads frc-ibl.2) and the PI model (frc-ibl.4) are
what replace them.

## Not included
Impurity or line radiation, ionization state (N2 rows assume Z=1, which is
invalid above roughly 10–20 eV), end and particle losses, rotation and the n=2
mode, plasma loading of the circuit, and switch or crowbar losses.

## Verification
`examples/frc_design_window.py` checks the closed forms independently:
- quadrature of the tanh profile against the flux formula, including zero net
  flux at the separatrix;
- a finite-difference null current against the decay formula;
- a sum of the repository's elliptic-integral loops against the sheet-coil
  centre field and the Wheeler inductance (flux linkage).
