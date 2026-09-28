# Initial Stage 0/1 delivery checks

This records the initial delivery. Subsequent dynamic MHD development and the
still-active overall objective are documented in `GOAL_STATUS.md` and
`MHD_VERIFICATION.md`; the initial 33-test count below is historical.

- 33 pytest checks passed on Python 3.13.7, macOS ARM64.
- The 18 vacuum/waveform checks passed before Stage 1 implementation began.
- CLI installed in the project virtual environment; help and package dependency checks passed.
- 3 normalized waveform experiments and a 27-case conductivity/amplitude sweep ran.
- Analytic/manufactured convergence, time-step, mesh and boundary-distance studies ran;
  measured errors are in VERIFICATION.md, with machine-readable values in outputs/verification/.
- Static field maps, comparison and phase-map images were visually inspected.
- Both HTML viewers embed 41 image frames, time sliders and playback controls,
  with external resource references removed. The 41-frame GIF decodes successfully.
- Live HTML browser interaction was not tested: the in-app browser's URL policy
  rejected local file navigation. File contents and embedded assets were checked
  without bypassing that restriction. Open the HTML files in a normal local browser.
- The repository was initialized with Git. No remote was added, no code was
  published and no commit was created. Large generated archives/viewers are ignored
  by Git; reproduction commands are in README.md.
- Every Python code file has a companion Markdown contract. Physics assumptions,
  SI conventions, source references and open research questions are documented.

Remaining limitations are scientific as well as numerical: fixed conductivity,
fixed plasma geometry and temperature, no momentum/energy dynamics, no Hall or
neutral physics, finite exterior boundary and unconverged topology durations.
The illustrative argon example strongly flags the missing-physics assumptions.
Stage 2 is a deliberate unimplemented interface, not a running MHD solver.
