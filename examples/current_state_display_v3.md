# Current signed-field presentation

## Purpose
Read-only successor to v2. Packages the four archived display frames with signed
current component controls, full-B field lines, current research status and a
scientific PNG containing quantitative signed meridional cross-sections.

## Components
- `research_inputs`: verifies the species checkpoint and its frozen sources,
  F16/F32 complete probe replays and F16 physical-audit replay, including raw
  files and non-runtime summary hashes. It does not run new scientific work.
- `package_currents`: packs `|J|, Jr, Jtheta, Jz` as float32 texture channels.
  Magnitude must reproduce the preceding texture bitwise. Signed components
  must round-trip exactly as float32; conversion errors against float64 are
  measured for every channel/time. Magnetic paths remain unchanged.
- `build`: assembles offline HTML, renders PNG, saves packing metrics and a
  manifest of every checked input. Refuses to overwrite an output directory.
- `main`: optional localhost-only serving through the preserved PNG handler.

## Contracts
No simulation, source repair, normalization of scientific states or time
interpolation. Grid and last evolved step remain F8/64 fluid, B32/64, step17 at
28.3333 µs. New diagnostic grids never replace these fields. The full physical
refinement audit is pending. Frame extrema are display samples, not certified
global extrema. Texture interpolation and illustrative optical blending are
distinct from the PNG's quantitative slice colors. Both include explicit units.

Full fresh-process builds must match HTML, PNG, packing metrics and manifest.
Static visual inspection and JS syntax checks do not establish browser UI
correctness when the Mac is locked. No unit tests. Freeze used sources; any
later presentation correction gets a separate named variant.
