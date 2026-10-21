# lowdensity_loss_check.py

A diagnosis of the 3D hybrid's low-density heat sinks (frc-ibl.11.9). It is
post-processing of saved snapshots only, with no state change. See the module
docstring for the method.

## Result on the v4 baseline shot (`outputs/hybrid3d_full_v4_uniform_h2`, 20 × 20 × 40 mm, 8 ions/cell, 5% floor)

**1. The sink sits where there are no carriers.**
- Below-floor ("thin") nodes hold 0.06–0.21 macro-ions on average.
- 100% of their removed heat comes from nodes with fewer than one macro-ion or
  with v_d > v_te.
- Nodes with no ions at all ("empty", inside the gas region) carry Joule
  heating from currents on their plasma edges. At those instants the solver's
  Spitzer resistivity is evaluated at the 1 eV maximum, so empty regions act as
  a resistive conductor, which the project rules exclude.
  - At 1 µs: empty 307 MW and thin 41 MW of Joule heat, against 26 MW in the
    resolved plasma.

**2. The size of the sink, from the solver's own counters.** By 4 µs it has
deleted 228 J by the cutoff and 144 J of electron energy on emptied nodes,
against 342 J of total Ohmic heat. Over 0–40 µs the removal is 43–60% of the
load work.

**3. Snapshot reconstruction overcounts.** Integrating snapshots at 0.1 µs gives
577 J removed over 0–4 µs, about 2.5× the counter. Snapshots catch empty nodes
just after their energy was deleted, at the 1 eV η maximum; between steps they
briefly hold transported heat. Use the counters for totals and the snapshots
for location.

**4. Physical counterfactual (thin nodes, J held fixed).**
- If the heat stayed in the electrons, Te would settle at 15–50 eV (median).
- The steady Spitzer dissipation, balanced by flux-limited parallel conduction
  (f = 0.1, axial connection length), is 10–40× below the numerical removal:
  e.g. 0.9 MW against 36 MW at 1 µs, and 3 J against 77 J integrated over
  0–4 µs.
- Removing the heat keeps Te cold and η maximal. This is a positive-feedback
  dissipation loop, not a loss model.

**5. Low-density amplification.** Below-floor nodes see 0 or 1 macro-ion per
step, so n, u = flux/n, the compression and the 1 eV-versus-hot η toggle are
PIC shot noise. This fits the CPU/GPU divergence that Codex measured: 13.7%
below the cutoff against 0.2% above it.

## Diagnosis
The removed heat is not standing in for physical losses. It is magnetic energy
dissipated by a numerical resistive medium in regions without resolved plasma:
the density floor in Ohm's law, plus the Spitzer η of a 1 eV, carrier-free node.
Physically:
- J would vanish there (no carriers), and the field would stay a vacuum field,
  storing energy rather than dissipating it.
- Where a thin plasma does exist, its electrons would heat, η would drop as
  Te^−1.5, and the heat would leave along open lines at ≲ 1/10 of the removed
  rate.

## Bounded experiment (run: see `vacuum_mask_replay.md`; it did not succeed, see the end)
A **dynamic vacuum mask** startup replay over 0–1 µs:
- Nodes with fewer than about one resolved macro-ion (or n < floor) are treated
  as exact vacuum for the field in each ion step: J = 0 on their edges.
- B there comes from the existing gap projector, applied to the moving mask.
- No floor carriers, no Joule heat and no removal in those nodes.

Compare against the v4 baseline over the same window:
- magnetic energy retained;
- Ohmic heat and removed heat;
- load work;
- the n = 0 FRC state;
- two-seed (or CPU/GPU) sensitivity in the formerly thin nodes.

**Success** means the removed heat falls to a few percent of the load work, the
retained magnetic energy rises by about the same amount, and the run-to-run
divergence below the cutoff drops to the above-cutoff level.

**Outcome.** Run as a `DeviceHybrid` subclass, the mask cuts the Joule-channel
sinks by 74%. Deleting the current of emptying cells costs more, 94 J by 1 µs,
so total removal rises from 68 to 112 J and retained magnetic energy falls by
41 J. The sink is the energy of current left behind by the departing ions. The
next test is ion-count convergence (ppc 8 → 32), not a different removal rule.

