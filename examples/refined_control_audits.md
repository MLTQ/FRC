# Independent audits of fully replayed refined controls

`run(destination, campaign_directory)` requires the complete serial physical
campaign, unchanged source/control records and complete same-source simulation
replays. For each accepted axial/radial small/large/tighter-bound control it reconstructs
every physical fluid/magnetic, electric and current component against the
existing same-grid epsilon1e-4 baseline with component_mode_audit. Every original
physical/reference fingerprint must match; each audit has a fresh full replay.
Strictly rejected controls are explicitly skipped, never assigned physical rates.

Each accepted small-amplitude comparison also gets a full independently replayed
component_axial_spectrum. Every one of its120 direct physical squared integrals
must exactly reproduce the separate component audit. Parseval and the complete
cell-mode spectrum remain; no field or state is filtered. These cell translation
modes fold continuous aliases and are not a complete unaliased continuum spectrum.

For each grid, all accepted new amplitude/residual controls receive current_rate_curl_audit and a
full replay. Independent strong curl acts inside polynomial cells and excludes
face distributions. Component errors and actual amplitude differences remain
distinct; total-current norms alone cannot certify weaker components.

Only one audit runs at a time, after the physical campaign. The parent records
active stage, source/campaign hashes, child audit paths, peak RSS and minimum
free disk observed at stage boundaries. Disk preflight is256MiB; these are not
continuous resource guards or separate accepted-solver tolerances. Each child
archive, exact raw integral array and partial stage survives. Fresh paths only.

This driver adds orchestration and cross-diagnostic comparison, not a new plasma
operator or independent physical reference. Individual full audit replays are
the verification evidence. No unit tests, numerical diffusion, filters, floors,
atmosphere, acceptance relaxation, time integration or formation claim.
