# hall_formation_compare.py

Compares a Hall (2.5D) formation run with a poloidal-only reference of the
same experiment, frame by frame at matching saved times.

```bash
python -m examples.hall_formation_compare HALL_RUN REFERENCE_RUN [OUT.json]
```

It reports, for each run:
- **End reconnection.** The first time from which the FRC stays formed for
  3 µs (closed, and shorter than 80% of the tube, as in the design scan). Also
  the first times the closed length falls below 80%, 60% and 40% of the tube.
- **Contraction.** The time and value of the minimum separatrix length.
- **FRC state.** Midplane r_s, l_s, trapped flux, ⟨n⟩ and Te+Ti at 2–20 µs.
- **Toroidal quantities.** The Hall run's peak B_θ, v_θ, J_r and J_z, and the
  B_θ and v_θ energies, from its history rows.
- **Cost.** Steps and wall time at the common end time. The step ratio is
  independent of machine load; the wall time is not.

It reads only saved outputs and evolves nothing.
