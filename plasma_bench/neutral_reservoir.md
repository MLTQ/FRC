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

## Optional charge exchange, escape limit (`cx_cross_section`)
- **Rate:** ions exchange with the reservoir at n_i K_cx, with
  K_cx = σ_cx √(8 e T/(π m_i) + u²) for neutrals at rest. σ_cx is declared;
  3·10⁻¹⁹ m² is the size of resonant H⁺–H exchange at 10–100 eV.
- **Shared neutral loss:** ionization and exchange draw on the same reservoir,
  consumed = n_n (1 − exp(−n_i (K_iz + K_cx) dt)), split by rate. At most half
  a cell's ions exchange in one step.
- **Each exchange:** a plasma ion is replaced by a cold ion at rest. The fast
  neutral leaves the system with that ion's bulk momentum, its bulk kinetic
  energy and 1.5 T of ion heat. Plasma mass is unchanged, and the neutral it
  came from is gone.
- **Why it is an upper bound:**
  - in a dense plasma the fast neutral's mean free path is millimetres, so most
    would be re-ionized or exchanged again inside the plasma;
  - escape is only right at the low-density edge;
  - for H2 targets at low energy, exchange is weaker than for H atoms.

## Bookkeeping
- `floor_mass` gains the ionized mass, and `floor_energy` loses the cost and the
  exchange loss. The
  solver's mass and energy residuals already subtract those closure sources, so
  they stay closed.
- Separate counters keep the ionized ions and mass, the total cost, and the
  loading heat.

## Not modelled (declared)
- Neutral motion, neutral pressure and wall interaction.
- Charge exchange beyond the escape-limit option: fast-neutral transport and
  re-ionization.
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

Charge exchange alone, at constant K_cx with n_i fixed:
- neutrals follow n_n0 exp(−n_i K t) to 1e-15;
- momentum follows p0 exp(−(n_n0/n_i)(1 − exp(−n_i K t))) at first order
  (observed order 1.00);
- mass is unchanged;
- the energy change equals the booked exchange loss to 1e-15.

