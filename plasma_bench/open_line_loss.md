# open_line_loss.py

## Intent
A **declared** energy-loss closure for field lines outside the separatrix.
Without it the gap plasma has nowhere to put its heat. In NUCTE-III formation
runs, Ohmic and numerical heating drove it to multi-keV temperatures and wall
pressures far above B²/2μ0: hot plasma held in by the glass rather than an FRC
held in by field. In a real theta pinch the edge is cooled by electron
conduction along open lines to the end walls.

## Model (Lie split after each accepted MHD step)
- **Which cells:** those with cell-averaged total flux ψ ≥ 0, i.e. open lines
  outside the reversed (ψ < 0) region.
- **Loss rate:** their thermal energy above the wall level relaxes on
  τ = 1.5 L/(f v_te), where:
  - L is the axial distance to the nearer tube end plus half a cell;
  - v_te = (e T_e/m_e)^½, with T_e = T_i = p/(2 n e);
  - f is `open_line_flux_limit` (~0.1, the standard free-streaming limit).

  This is the flux-limited conduction loss of a column of length L.
- **Integration:** exact, with τ frozen over the step.
- **Wall level:** n·2e·T_wall/(γ−1), with T_wall = `open_line_wall_temperature_ev`.
- **Cold dual-energy cells** update ρA consistently.
- **Accounting:** removed energy is reported as `open_line_loss_j` and added
  back in both energy residuals.

## Limits
- It is a sink model, not anisotropic conduction: no transport along lines
  inside the domain, no conduction across lines, and no heating of the ends.
- Before reconnection the reversed bias flux (ψ < 0) counts as closed, even
  where it still reaches the ends.

## Verification (`python -m examples.open_ends_verification --conduction`)
In a uniform resting plasma where every line is open, the solver matches an
independent ODE integration of dU/dt = −(U − U_w) f v_te/(1.5 L):
- mid-domain cells agree to 6e-5 and 1e-3;
- the end-adjacent cell is off by 11%, because the cooling gradient drives
  flow and compression there;
- the energy residual is 2e-10 J.

2.5D: the v_θ energy joins the kinetic energy. The caller's `magnetic`
includes the B_θ energy.
