# Latest fields with quantitative cross-sections

Presentation-only successor to the frozen v2 scientific PNG. `render` receives
verified float64 display arrays and emits the last accepted state with two 3D
panels: all-component current magnitude and full magnetic field lines. It adds
unshaded meridional `Jz` and `Bz - 5 mT` images with signed physical legends,
making structure visible without increasing magnetic line deflection.

All plotted values derive directly from the same saved frame. The 3D geometry
preserves physical aspect ratio and includes the full bias field. Cross-sections
use cell-center display samples, constant edge extension over the outer half
pixel, and nearest-pixel rendering. These are display approximations, not new
independent physical samples or solver convergence evidence. Current sign uses
the physical +z axis. Axial periodicity and the filled domain remain explicit.

Labels preserve accepted step 17, 28.3333 µs, rejected step 18 and the absence of
demonstrated FRC formation. The current research caption distinguishes scalar
source refinement from unresolved full-system growth. Full build replay and
visual PNG inspection are required; no unit tests or simulation are run.
