# nucte_compare.py

`python -m examples.nucte_compare RUN [OUT]` reads a `nucte_formation` run's
`snapshots.npz` and `history.json`. At every 1 µs snapshot it computes:
- the separatrix profile r_s(z): the zero of the total flux off a reversed
  axis, infinite where the reversal reaches the wall;
- the reversed segment containing the midplane, and whether it is closed
  (unreversed or wall-free boundaries at both ends);
- the trapped poloidal flux, 2π|ψ_min|;
- the volume-averaged density and Te+Ti = p/(n e) inside r < r_s(z), assuming
  Z = 1 and one temperature;
- the midplane side-on chord integral ∫n dl.

It writes `comparison.json` with the NUCTE-III case-1 targets (Takahashi 2000):
- r_s(0) = 6.2 cm, l_s = 62 cm;
- axial contraction 6.5 µs at 3.8 cm/µs;
- lifetime 65 µs, n=2 onset at 28 µs;
- n̄ = 2.2e21 m⁻³, Te+Ti = 120 eV.

Those are equilibrium values with unpublished shot waveforms, so comparisons
are indicative.

`python -m examples.nucte_figure RUN [RUN2 ...] OUT.png` draws:
- r_s(z) at 2–15 µs for the first run;
- small multiples of the midplane radius, length, density, temperature and
  trapped flux, each against the published case-1 value.
