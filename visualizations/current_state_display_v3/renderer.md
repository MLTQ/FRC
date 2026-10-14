# Signed current and full magnetic renderer

Presentation successor to the frozen trajectory renderer. The same camera,
geometry, field paths and volume quadrature are retained. Texture RGBA now means
`|J|, Jr, Jtheta, Jz` in A/m². `state.component` selects the displayed channel.
Magnitude uses the warm sequential ramp; signed components use blue / neutral /
orange. Opacity depends on absolute component strength. This is illustrative
optical integration, so composite pixels are not quantitative surface samples.

Shared or per-frame limits come from the archived float64 arrays. Signed limits
are symmetric about zero. No smoothing, field deformation, physical evolution,
current streamlines or changes to magnetic paths are introduced. Linear texture
interpolation, radial edge clamping and cutaway limitations remain those of the
original trajectory display. Browser interaction must be checked separately;
static build/replay does not verify WebGL behavior. Preserve once used.
