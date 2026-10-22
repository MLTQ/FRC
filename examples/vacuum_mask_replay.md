# vacuum_mask_replay.py

Bounded experiment for bead frc-ibl.11.9. It asks whether the 3D hybrid's
low-density heat sinks (`examples/lowdensity_loss_check.md`) are set by a
numerical resistive medium rather than by physics. The method is to replay the
v4 baseline startup with carrier-free cells treated as exact vacuum
(`plasma_bench/hybrid3d_vacuum_mask.md`).

## Runs
All runs use the v4 baseline configuration (`outputs/fine_uniform_h2` from
t = 0):
- 20 × 20 × 40 mm cells;
- 8 ions per cell;
- 5% density floor;
- coupled 3D banks;
- conservative energy electrons with Spitzer η;
- dt = 1 ns, to 1 µs.

They go through `examples.hybrid3d_frc.main` unchanged. Only the model class
(and extra history fields) is swapped.

| run | seed | mask |
|---|---|---|
| `outputs/hybrid3d_full_v4_uniform_h2` (existing, cut to 1 µs) | 0 | no |
| `outputs/vacuum_mask_v4_base_seed1` | 1 | no |
| `outputs/vacuum_mask_v4_seed0` | 0 | yes, min_ions = 1 |
| `outputs/vacuum_mask_v4_seed1` | 1 | yes, min_ions = 1 |

`compare` writes `outputs/vacuum_mask_v4_seed0/vacuum_mask_compare.json`.

## Check
While no cell has emptied, the masked seed-0 run is bit-identical to the
baseline: node density at 0.1 µs differs by exactly 0.

## Result at 1 µs (load work 771–774 J in all four runs)
| | baseline s0 | baseline s1 | mask s0 | mask s1 |
|---|---|---|---|---|
| magnetic energy (J) | 540.4 | 539.0 | 499.6 | 498.0 |
| response magnetic (J) | 384.0 | 382.3 | 338.6 | 336.9 |
| Ohmic heat (J) | 95.4 | 94.7 | 48.1 | 47.9 |
| cutoff-removed heat (J) | 47.8 | 47.1 | 6.0 | 6.0 |
| empty-node energy removed (J) | 20.2 | 20.4 | 11.5 | 11.6 |
| mask-change projection loss (J) | – | – | 94.5 | 93.6 |
| **all numerical sinks (J)** | **68.0** | **67.5** | **112.0** | **111.2** |
| load-energy residual (J) | −0.03 | 0.00 | −94.6 | −93.7 |
| masked gas cells (of 9900) | – | – | 3793 | 3790 |
| separatrix at midplane (m) | 0.121 | 0.120 | 0.114 | 0.113 |

Trapped flux is 1.6 mWb in all four runs and none is closed yet. The residual
of the masked runs equals the projection counter to within 0.05 J, so their
accounting closes.

**Seed-to-seed B_response difference** (relative RMS, at matched times):

| t (µs) | baseline plasma / thin / empty | mask plasma / thin / empty |
|---|---|---|
| 0.5 | 0.27 / 0.32 / 0.42 | 0.25 / 0.29 / 0.24 |
| 0.8 | 0.22 / 0.47 / 0.41 | 0.20 / 0.33 / 0.17 |
| 1.0 | 0.15 / 0.53 / 0.31 | 0.14 / 0.31 / 0.14 |

## Reading
- **The mask does not succeed by the stated criterion.**
  - The Joule-channel sinks fall from 68 J to 17.5 J (8.8% → 2.3% of load work).
  - But deleting the current in cells as they empty costs 94 J. The total
    numerical removal rises to 112 J (14.5%), and the retained magnetic energy
    falls by 41 J instead of rising.
- **Why.** By 1 µs the radial implosion has emptied the outer 38% of the gas
  cells. The trailing part of the current sheet sits in cells that are losing
  their last macro-ions. The baseline dissipates part of that carrier-free
  current resistively. The mask deletes all of it at once, as the
  minimum-energy projection does.
- **What this means for the diagnosis.**
  - It confirms that the low-density sink is energy held in current with no
    resolved carriers.
  - It shows that the sink is set upstream, by the field decoupling from the
    departing ions at the floor density. How that current is removed does not
    set it.
  - Physically, the sheet current would move inward with its carriers (J × B
    work on the plasma). It would not be left behind to be dissipated or
    deleted.
- **Noise.** The mask halves the seed sensitivity in thin and empty nodes
  (0.53 → 0.31 and 0.31 → 0.14). Thin nodes are still twice as noisy as the
  plasma (0.14). The divergence target is not met either.

## Ion-count convergence (no mask; `run ... mask=0 ppc=N`)
The loader places ions on a side³ lattice, so `ppc` rounds to a cube: 16 and 32
both give 27. Velocities are the only seeded quantity. The two 27-ion seed-0
runs (`outputs/ppc_conv_v4_p27_s0_repeat`, `outputs/ppc_conv_v4_p32_s0`) are
bit-identical.

All values at 1 µs, seed 0 (the seed-1 value in brackets):

| ions per cell | 8 | 27 | 64 |
|---|---|---|---|
| cutoff + empty-node removal (J) | 68.0 (67.5) | 53.3 (54.8) | 45.2 |
| Ohmic heat (J) | 95.4 (94.7) | 82.3 (83.4) | 75.9 |
| magnetic energy (J) | 540.4 (539.0) | 548.1 (548.4) | 550.9 |
| kinetic energy (J) | 263.2 | 261.0 | 260.6 |
| separatrix at midplane (m); trapped flux | 0.121; 1.6 mWb | 0.121; 1.6 mWb | 0.121; 1.6 mWb |
| macro-ions per thin node | 0.16 | 0.50 | 1.13 |
| empty gas nodes (snapshot) | 3788 (s1) | 3322 | 2912 |
| Joule power, empty / thin / plasma nodes (MW, snapshot) | 321 / 40 / 28 (s1) | 268 / 45 / 23 | 235 / 45 / 21 |
| seed-to-seed B difference, plasma / thin / empty | 0.15 / 0.53 / 0.31 | 0.12 / 0.51 / 0.28 | – |

Snapshot rows come from `lowdensity_loss_check`. The seed differences come
from `seed_difference`.

- **The removal falls with carriers, but slowly.**
  - The successive drops are 14.7 J and 8.1 J. A power law through three
    points gives about N^−0.25, half the shot-noise exponent of N^−0.5.
  - Its asymptote (about 10 J) is poorly constrained by three points.
  - At 64 ions per cell, 45 J (5.9% of the load work) is still removed by 1 µs.
- **The formation state at 1 µs is already insensitive to the ion count.**
  Separatrix and trapped flux agree to 1%, and the magnetic energy changes by 2%.
  The sink is a slowly converging energy-accounting error, not yet a change in
  the configuration.
- **Most of the sink is Joule heat in genuinely emptied nodes** behind the
  imploding sheet.
  - Their number falls only 23% from 8 to 64 ions per cell, and they carry
    80–85% of the snapshot Joule power.
  - The floor closure treats them as a 1 eV Spitzer resistive medium. Its
    magnetic diffusion time across a 20 mm cell is about 0.5 µs, comparable to
    the implosion, which is about the worst case for dissipation.
- **The seed divergence below the cutoff is not shot noise.** Thin-node
  divergence barely changes (0.53 → 0.51) while carriers per node triple. It is
  dynamical sensitivity of the trailing edge, which fits the CPU/GPU roundoff
  amplification Codex measured.

## Edge-model replay (`run ... mask=2 [min_ions=N]`; `plasma_bench/hybrid3d_edge.md`)
All values at 1 µs, seed 0 (the seed-1 value in brackets):

| | floor (base) | mask | edge, 1 ion | edge, 2 ions | edge, 4 ions |
|---|---|---|---|---|---|
| total numerical removal (J) | 68.0 (67.5) | 112.0 (111.3) | 103.8 (101.9) | 121.2 | 132.7 |
| mask-change projection (J) | – | 94.5 | 93.5 (91.6) | 107.8 | 116.5 |
| empty-node electron energy (J) | 20.2 | 11.5 | 4.7 (4.5) | 8.9 | 13.5 |
| cutoff-removed heat (J) | 47.8 | 6.0 | 5.7 | 4.5 | 2.7 |
| electron energy (J) | 33.4 | 39.9 | 42.3 (43.1) | 44.1 | 43.7 |
| magnetic energy (J) | 540.4 | 499.6 | 492.6 (491.9) | 487.7 | 486.2 |
| pressure-centering residual (J) | 0.03 | 0.04 | −5.0 (−5.2) | −5.2 | −4.2 |
| separatrix at midplane (m) | 0.121 | 0.114 | 0.113 | 0.113 | 0.113 |
| wall time (s) | 153 | 136 | 194 | 260 | 268 |

- **The electron edge works here too.** It cuts the edge electron-energy removal
  from 11.5 to 4.6 J.
- **The current deleted at emptying still dominates, and it grows with the
  vacuum threshold.** Here the sparse plasma (1–4 macro-ions per cell) carries
  the formation current. Declaring it vacuum deletes that current, while the
  floor closure dissipates it more slowly.
- **The edge model does not improve the startup.**

## Carrier ledger and electron edge without a mask (`run ... mask=0 ledger=1 electron_edge=0|1`)
`plasma_bench/hybrid3d_carrier_ledger.md`; all values at 1 µs.
- **Checks.** With `electron_edge=0` the run reproduces the original baseline:
  identical counters, n and p_e, and B within 9e-13 T. The same-seed repeat of
  the electron-edge run is bit-identical. The ledger class sums match the solver
  counters to about 1e-13 J.

| | baseline s0 | electron edge s0 (= repeat) | electron edge s1 |
|---|---|---|---|
| Ohmic heat, orphaned / sparse / resolved edges (J) | 43.1 / 2.9 / 49.4 | 68.4 / 2.2 / 45.8 | 68.0 / 2.2 / 45.7 |
| cutoff removal, thin / empty nodes (J) | 12.7 / 35.1 | 16.8 / 51.4 | 16.7 / 50.8 |
| empty-node electron energy removed (J) | 20.2 | 6.0 | 5.9 |
| **total numerical removal (J)** | **68.0** | **74.2** | **73.4** |
| magnetic energy (J) | 540.4 | 518.5 | 517.5 |
| electron energy (J) | 33.4 | 46.1 | 46.9 |
| \|J\| A on orphaned / sparse / resolved edges (kA) | 8.4 / 3.0 / 12.5 | 5.6 / 2.4 / 11.9 | 5.8 / 2.3 / 11.9 |
| pressure-centering residual; load residual (J) | 0.025; −0.025 | 0.039; −0.002 | 0.042; 0.022 |
| separatrix (m); trapped flux | 0.121; 1.6 mWb | 0.120; 1.6 mWb | 0.120; 1.6 mWb |

At 0.5 µs, orphaned edges carried only 0.4 J of the 27.5 J of Ohmic heat.
About 5,700 cells lost their last ion during 0–1 µs.

- **Where the left-behind current's energy goes** (floor closure):
  - By 1 µs, 35% of the current magnitude sits on orphaned edges, where no
    adjacent cell holds an ion. It dissipates 43 J there, 45% of all Ohmic heat.
  - The heat is then deleted: 35 J by the cutoff on empty nodes, plus the
    empty-node energy removal.
  - So the startup sink is mainly Ohmic dissipation of carrier-free current,
    which forms as the implosion empties the outer gas after 0.5 µs.
  - Resolved-thin current (the sparse class) contributes only 2–3 J.
- **The electron edge alone does not justify itself.**
  - It blocks the leak: empty-node removal falls from 20 to 6 J and the electrons
    keep 13 J more.
  - But Ohmic heat on orphaned edges rises from 43 to 68 J. Total removal rises
    from 68 to 74 J, and magnetic energy falls by 22 J.
  - Inferred mechanism: `ThermalElectrons.temperature` divides by
    max(n, 1 m⁻³). Electron energy that leaked into an empty node within a step
    gave it an enormous Te and almost zero Spitzer η until `finish_step`
    removed it. Blocking the leak exposes the 1 eV η cap on the orphaned
    current.
  - **The baseline's lower orphaned dissipation is therefore an artifact of the
    leak.** η at a node with no carriers is undefined, and the closure either
    zeroes it by accident or caps it at 1 eV.
- **Answer to the coordinated question.**
  - Electron-edge transport cannot be justified on its own. It is correct for
    the electrons, but it makes the field dissipation of the carrier-free current
    worse.
  - When carriers leave a cell, their current stays in the field. No ion is left
    to receive J × B work, so the current can only dissipate, at an η the
    closure does not define, or be deleted (the mask).
  - Any physical edge model has to stop carrier-free current forming, by keeping
    the field coupled to departing ions. Choosing how to remove it afterwards
    cannot fix this.

## Floor convergence with the carrier ledger (`run ... mask=0 ledger=1 electron_edge=0 floor=F`)
Base closure; all values at 1 µs, seed 0 (the seed-1 value in brackets).

| Ohm's-law floor (fraction of peak n) | 5% | 2% | 1% |
|---|---|---|---|
| **numerical removal, cutoff + empty-node (J)** | **68.0** | **48.4** | **39.5 (40.5)** |
| cutoff removal (J) | 47.8 | 22.1 | 11.8 (12.0) |
| empty-node electron energy (J) | 20.2 | 26.3 | 27.6 (28.6) |
| Ohmic, orphaned / sparse / resolved edges (J) | 43.1 / 2.9 / 49.4 | 21.0 / 3.5 / 52.4 | 12.2 / 4.2 / 57.0 |
| \|J\| A, orphaned / sparse / resolved (kA) | 8.4 / 3.0 / 12.5 | 7.6 / 3.5 / 13.4 | 6.6 / 3.3 / 12.9 |
| wall electron enthalpy outflow (J) | 38.7 | 55.5 | 73.6 (73.2) |
| magnetic / response magnetic energy (J) | 540.4 / 384.0 | 518.4 / 356.9 | 509.0 / 346.4 |
| ion kinetic / electron energy (J) | 263.2 / 33.4 | 272.7 / 35.7 | 271.2 / 34.6 |
| load work (J); load residual (J) | 776.0; −0.03 | 758.1; −0.10 | 754.1; −0.19 |
| **separatrix at midplane (m)** | **0.121** | **0.112** | **0.112 (0.111)** |
| trapped flux (mWb) | 1.6 | 1.6 | 1.5 (1.6) |
| field substeps; wall time (s) | 3508; 145 | 7510; 228 | 14186; 363 |

The ledger class sums match the solver counters to 5e-13 J. At 1%, seeds 0
and 1 agree within 2.5%.

- **The floor sets the dissipation of carrier-free current.**
  - Orphaned-edge Ohmic heat falls from 43 to 21 to 12 J as the floor goes from
    5% to 2% to 1%, while the orphaned current magnitude falls only 21%.
  - A lower floor makes the empty medium Hall-dominated (faster electron
    response), so its current dissipates less. This is the mechanism the
    previous replays pointed to.
- **Total numerical removal falls as about floor^0.65:** 68 → 48 → 40 J. A
  three-point extrapolation to floor → 0 gives about 24 J (3% of load work),
  poorly constrained.
  - What remains is mostly electron energy deleted on emptied nodes, which
    rises from 20 to 28 J as the floor drops.
  - The wall enthalpy outflow also rises (39 → 74 J), because the floored Hall
    velocity in thin edge nodes is larger. It is a declared boundary channel,
    but it is not converged either.
- **The formation state is not converged at 5%, unlike in ion count.** The
  separatrix at 1 µs moves from 0.121 to 0.112 m (−8%) going to 2%, then agrees
  at 1%. Trapped flux stays at 1.6 mWb.
- **Design consequence.** The 5% production floor biases the early separatrix
  radius by about 8%. A 2% floor converges it at 1 µs for about 2× the field
  substeps (1.6× wall time). Whether that holds beyond 1 µs is untested.
- **Next bounded step (not run):** 1% floor plus the electron-edge block. At a
  low floor the orphaned dissipation that the block exposes is 3.5× smaller, so
  the block may now pay off against the 28 J of empty-node removal.

## Where this leaves frc-ibl.11.9
The sink comes from the sheet current left in emptied gas, which the floor
closure turns into a resistive conductor.
- More carriers shrink it slowly.
- Deleting it at emptying (the mask) costs more than dissipating it.
- Making the emptied region more conducting would violate the exact-vacuum rule.

The sharp-edge model (above) confirms this from the other side. Removing the
floor medium deletes the current instead. Resolving the current-carrying sparse
plasma is the remaining route. Until then:
- report the removal as an accounting uncertainty that converges with ions per
  cell (about 6–9% of load work at 1 µs for 64–8 ions per cell);
- treat the formation geometry at 1 µs as converged in ion count.
