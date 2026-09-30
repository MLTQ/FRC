# Local scientific volume and field-line renderer

WebGL2 ray integration displays the archived |J| texture revolved about its
axis. Color and optical opacity are presentation choices; there is no optical
emission model. The removable quarter-sector is a display cutaway only. Float
texture interpolation is explicit and fixed color limits span all snapshots.

Thin six-sided tubes follow measured poloidal flux contours in meridional
planes. Toroidal rings indicate the actual current direction geometry; no
helical/toroidal magnetic component is invented. Width, glow and density of
seeded lines have no quantitative field-strength meaning. Their color does.
The camera uses perspective and a synchronized orbit for comparison panels.

`FieldView` owns one canvas/context and replaces its buffers between saved
frames. `draw` supports higher-resolution export. No network or third-party
runtime is required. The accompanying page handles controls, labels and export.

The default camera fits both the volume and exterior lines below the panel
labels. Current tubes obey the same cutaway as the volume. Five meridional
planes show closed lines; three show exterior lines, avoiding an opaque cage.
