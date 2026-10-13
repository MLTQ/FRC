# Additional camera controls

Extends the frozen trajectory controller after it initializes shared `state`
and `dirty`. Presets provide perspective, transverse side and axial end views
using its existing camera convention (physical z is scene x). Each preset
stops camera orbit. Keyboard arrows rotate; plus/minus zoom within the original
camera limits. Both canvases remain synchronized and become keyboard focusable.

No new field sampling, time interpolation, deformation or physics. These
controls need browser end-to-end interaction to claim UI verification; static
generation and JavaScript parsing alone do not establish that.
