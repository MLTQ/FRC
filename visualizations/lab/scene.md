# scene.js

## Purpose
Dependency-free perspective 3D renderer of saved physical electron/ion positions,
full cage edges and chamber reference circles. Mouse/touch orbit and zoom change
only the camera. No rotational particle copies or fabricated trails.

## Components and contracts
- `Scene`: geometry/frame inputs, display-only count, camera reset, draw. A
  stratified species subset keeps electrons and ions visible at reduced counts.
  Points are depth sorted; steel wire line width is illustrative, not a mesh.
- `drawSlice`: actual saved y=0 Poisson cut, signed physical voltage scale, null
  outside gas. Cell coloring interpolates no new solver resolution.
- `drawHistory`: actual kinetic and electric field energy in nJ, saved timeline.
- `drawDrive`: saved live-ion RMS radius (left, mm) and actual cage voltage
  (right, kV) against microseconds, with a selected-state marker. No curve is
  invented before two states; an absent cloud breaks the radius trace.
  Radius/voltage use distinct labelled scales; no predicted frequency appears.
- Microsecond experiment timelines use µs; older short runs retain ns.
- Optional `paths` and `samples` render saved E/diagnostic B lines and current
  glyphs in the same freely orbiting 3D geometry. Lines carry direction arrows;
  colors map log magnitude, and current glyph lengths are display-normalized.
  `fields.js` supplies source/scope/physical scales. `canvasContext` and `color`
  are shared by the new physical moment and phase plots.

CPU Canvas projection is a true interactive 3D view; it is not a GPU solver.
Frame sampling changes dot visibility and provides no particle convergence proof.
Canvas backing stores resize only when their physical pixel dimensions change;
each draw resets its transform without reallocating/blanking the backing store.
This also avoids resize-observer churn during playback.
