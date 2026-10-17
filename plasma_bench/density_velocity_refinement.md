# Mechanical nodes at represented density jumps

`refine_density_jumps` adds existing fine density faces to the mechanical mesh
when adjacent positive cell densities differ by at least an explicitly supplied
factor. It preserves every inherited velocity knot, every fine density cell,
support and population. The criterion uses differences of log(N_relative/h),
so a common physical number scale cancels without forming underflowing densities.
Geometry differences are evaluated with the selected precise backend first.

The contrast is a disclosed spatial-resolution control, not a physical rate or
a population cutoff. No cell is dropped or assigned a different density, heat,
mass or force. This helper returns knots and diagnostics only; callers must
represent the existing velocity field or re-solve the declared conservative
source transfer on that richer space. Actual new acoustic bounds still apply.
There is no spacing relaxation, close-node merge or automatic density cap.

The motivating complete neutral Euler experiment keeps65 fine density cells
with a1e8 density jump inside one of64 velocity elements. It produces a3040.63eV
peak and fails when a334km/s tail reaches the remote boundary at0.638us, despite
energy conservation near1e-11. Adding just the represented jump face gives65
velocity elements, completes the original1.151us and has peakT0.9253eV. The
independent exact shock maximum is1.2042eV; accuracy still requires refinement.
See outputs/material_density_step_interface_comparison for complete histories.

The completed actual-stage trial now supports correction of the Ar9 spike.
The explicit contrast1e4 trial adds4 knots (38→42 elements) to its captured
post-source initial state; every original N/U and physical velocity field is
retained. The smallest mechanical spacing is1.34581um. Full5ns transport
completes in221 steps with Erel5.10e-12; final Ar9 maximum falls370.16→0.46843eV
and volume RMS excess over reversible compression falls29.73→0.01573eV.
Ar8 remains14.176eV. Independent80-digit saved-coordinate readout and inspected
plots are in reaction_density_jump_transport. Reversible compression is not an
exact shock reference; convergence and other heavy-species errors remain open.
A separate frozen full20ns coupled candidate is still running, and a preserved
all-heavy stage remains queued. No production source default is changed.

Inputs are ordered precise density faces, strictly positive finite relative
amounts, existing ordered nested velocity knots and a finite contrast>1. Output
is their ordered union with selected faces plus selection/count/contrast records.
