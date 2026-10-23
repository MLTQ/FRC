# neutral_reservoir.py

A static neutral reservoir for formation from a partially pre-ionized fill
(bead frc-ibl.4). `NeutralReservoirMHD` subclasses `ConservativeVacuumMHD`
without changing it. After each accepted MHD step (and after the density floor
and open-line loss), it applies `apply_ionization` as a Lie split, in the same
way `density_floor.apply_floor` works.

## Equations, per fluid cell, over one step dt
- **Ionized atoms:** dN = n_n (1 − exp(−n_e K_iz(T) dt)). This is exact for
  n_e and T frozen over the step. n_e = n_i = ρ/m_i. T is the single-fluid
  temperature from p = 2 n T, using the cold-cell entropy energy where dual
  energy marks a cell cold.
- **Rate:** K_iz(T) is the NIST-BEB Maxwellian total ionization of ground-state
  H2 (`molecular_beb`), interpolated in log–log over 0.3–2000 eV and zero below
  0.3 eV.
  - It is applied per atom-equivalent of the reservoir. Each event is assumed to
    become two protons promptly, because H2⁺ dissociation and atomic ionization
    are taken as fast at the sheath temperatures that matter.
  - This is declared, not computed: neither H atoms nor D2 are modelled.
- **Mass and momentum:** the new ions are born at rest, so mass rises by m_i dN
  and momentum is unchanged. The mass-loading kinetic-energy drop stays in the
  cell as heat, so total energy is unchanged by it.
- **Cost:** the declared ionization cost E_c per new ion (dissociation share,
  ionization, excitation and radiation) is removed from thermal energy. It is
  capped so a step never spends more than half a cell's thermal energy.

## Bookkeeping
- `floor_mass` gains the ionized mass and `floor_energy` loses the cost. The
  solver's mass and energy residuals already subtract those closure sources, so
  they stay closed.
- Separate counters keep the ionized ions and mass, the total cost, and the
  loading heat.

## Not modelled (declared)
- Neutral motion, neutral pressure and wall interaction.
- Charge-exchange drag and charge-exchange loss of fast neutrals.
- Electron–neutral resistivity.
- Recombination.
- Molecular and atomic stages, and a computed E_c.
- The ionization source only adds plasma where electrons already are. It never
  creates plasma in exact vacuum: the exterior is not part of the fluid grid.

## Verification
`python -m examples.partial_pi_formation check` compares the operator with the
exact logistic solution of dn_e/dt = K n_e n_n at constant K. Over 3 e-folds:
- convergence is first order (0.98, 0.99, 0.99);
- mass is conserved to 1e-16;
- energy changes by the cost to 1e-15;
- a moving cell keeps its momentum exactly, and its energy changes only by the
  cost.
