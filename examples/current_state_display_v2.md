# Latest field viewer with the corrected scientific poster

## Purpose
Presentation-only successor to `current_state_display.py`. Uses
`current_state_poster_v2.render` to correct clipped geometry in the first PNG.
The interactive viewer is identical; no field, camera-control or physics
equation changes. The initial output directory is retained unchanged.

## Components and contracts
`verify_inputs`, `replace_once`, `build` and `main` keep the contracts of
`current_state_display.md`. The default destination is
`outputs/current_state_display_v2`; completed output cannot be overwritten.
The manifest records the corrected renderer and its companion. The localhost
server, when requested, uses port 8787 and serves only the selected directory.

All four original saved frames and audit evidence must retain their hashes.
The last accepted time remains 28.3333 microseconds. New audits do not create
a new physical state. Compare a full fresh-process replay and visually inspect
the PNG; browser interaction remains separate from generation verification.
No unit tests. Preserve this source after executing it.
