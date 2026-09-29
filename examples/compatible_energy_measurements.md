# Compatible magnetic-energy measurement protocol

Runs the unchanged independent shock through nr=384, source pulses and analytic
smooth references using staggered energy, natural vacuum traces and SSP RK3.
RK2, half-CFL, longer chamber, isolated energy/trace ablations, disabled limiting
and farther magnetic matching surfaces expose method and boundary sensitivity.
This is a fixed-chamber verification study, not a formation experiment.

Cases select every new method control explicitly. Package/common measurement,
driver and environment fingerprints guard reuse; input options must match too.
Each final record is persisted separately. During expensive runs, accepted-state
summaries are printed and saved with a `.progress.json` suffix. These records
say `completed=False, status=running` and cannot be reused as final results.
They provide observability, not numerical restart data. The final result removes
its temporary progress file. Failed cases retain their last accepted state.

`--case` selects independent named cases without concurrent writes to one case.
The historical cell/Riemann/RK2 studies and their artifacts remain separate.
The reporting module owns plots and interpretation; edits there do not invalidate
these numerical measurements. Accepted Euler-stage counts, rather than twice
the step count, normalize RK2/RK3 limiter activity.
