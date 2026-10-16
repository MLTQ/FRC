# Retained-resolution moving argon and actual cost

Requests the unchanged full 20 ns neutral/ion/electron argon experiment with
physical masses, 1e14/m³ per initially present component, 10 mm half-width,
0.1/0.2/2 eV initial temperatures and the actual ground-ladder DI/RR source.
Array-assembled moments change computational assembly only. Source-half /
transport-full / source-half and half the current material timestep are retained.

`minimum_density_width=0` forbids physical-width density coarsening. The existing
explicit 8-ULP cleanup remains, with all species support ends protected. This
is not a claim that ULP cleanup leaves local current exactly invariant. Requested
mechanical resolution defaults to32; `--cells` changes it with the original
plasma-phase scaling. No density-cell-count limit is used.

Measures each source half and transport step separately, actual density-cell
growth, peak process memory and disk headroom. The declared 600 s observation
budget and 2 GiB peak-memory bound stop only between coupled steps. They are
run budgets, not accuracy criteria or proof of exhausted machine resources.
Each accepted state is written immediately using existing normalized N/U records;
the summary carries actual time, populations, photons and the complete physical
energy. A conservative uncompressed-size estimate reserves256 MiB free storage.
If the accepted state cannot be saved within that reserve, this is recorded
explicitly; its field history is not invented. There is no shortened success.

Source or geometry rejection preserves prior accepted snapshots and the actual
failure message. No retries, floors, species removal, mass changes or heat repair.
The scalar source comparison and independent complete conversion histories are
in examples/atomic_vector_moments.md. This run does not prove resolved current,
physical closure validity or any formation requirement by itself.

Use `python -m examples.reactive_material_fine FRESH_OUTPUT [--cells 32]`.
