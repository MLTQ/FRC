# Detailed volume display from completed saved states

`build` augments an existing experiment display without retracing fields or
running the solver. It preserves every actual timestamp and links its original
frame files and scientific figures. A fresh destination is mandatory.

`volume` evaluates the actual compatible material basis on 65×257 strictly
interior labels, using the experiment's declared initial density. J=enw is
mapped to a 96×384 regular physical (r,z) texture by linear interpolation.
The fourth channel records interpolation support. Empty exterior remains empty;
the renderer only samples cells with four supported corners. This conservatively
omits a thin boundary strip. The texture is a display approximation, not a new
solver mesh. Spatial derivative jumps can be visually interpolated.

The full physical response supplies drift coefficients. Interior sampling uses
only its space/fluid view, avoiding unused prescribed-coil B evaluation at every
volume label. Current has no additive coil-source term; all coil effects in the
coupled response are retained. No physical-model files are changed.

Each little-endian float32 texture contains Jr,Jtheta,Jz in A/m² and support.
Texture bounds are physical (r,z) divided by the original display extent.
Per-state and shared maxima include dense basis samples and the earlier cloud;
there is no renormalization of the underlying data. Input hashes are retained.

Run `python -m examples.material_volume_display EXISTING_DISPLAY FRESH_OUTPUT`.
Generated textures and linked original frames remain local and Git-ignored.
`--limit N` is for an explicit partial display build during UI development;
it does not claim a complete playback. The final deliverable uses all frames.

The builder injects volume.js before the shared viewer.js into volume.html.
Both views share a physical orthographic camera. Field paths may be rotated
azimuthally because the underlying solution is axisymmetric, never distorted.
No new physical 3-D modes, temporal interpolation or instability are introduced.
