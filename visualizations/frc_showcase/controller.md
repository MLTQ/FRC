# Viewer interaction and labeled poster export

Controls operate only on three archived snapshots and presentation camera state.
Dragging either panel rotates both; wheel zoom is bounded. Time selection
replaces measured arrays without interpolating synthetic intermediate times.
Auto orbit animates the camera only. Cutaway removes one viewing quadrant.

Poster export re-renders both panels at high resolution, adds scientific units,
fixed scales, saved time, physical size and interpretation. Served on localhost,
the PNG is posted to the narrow `/save-poster` endpoint. A standalone file uses
a normal browser download instead. The page has no other network activity.

The initial view uses a quarter cutaway and oblique camera to reveal current
structure. Reset restores the same fitted camera; cutaway stays user-controlled.
