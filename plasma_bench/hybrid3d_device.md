# hybrid3d_device.py

For known coil sources, plasma current is curl(B−B_source)/μ0. The known
source curl is subtracted at each RK stage before Ohm's law, heating and
particle forces. Otherwise an unresolved coil stencil can spuriously push
cold gas. This subtraction retains the full total B in the Lorentz field;
it removes no plasma mode or response current.

Plasma current uses the same edge support as plasma E. Conductor sheet
currents on excluded edges cannot enter the node Hall force. Edge current
averaged to a node is divided by its gas dual-volume fraction, matching the
normalization of PIC density and flux. With these supports, Hall J·E is
zero and ideal-field work has the same volume weight as ion acceleration.
Previously unmatched conductor currents injected artificial magnetic energy
in full-shot startup; those trajectories are invalid.

The 3D hybrid model in the 12" device: milestone M2 of
[docs/HYBRID_3D.md](../docs/HYBRID_3D.md). It shares the particle kernels and
Ohm's law of [hybrid3d.py](hybrid3d.md). The field kernels are numba versions,
identical to the numpy operators to round-off.

The device field-step bound uses the exact Yee difference-symbol norm
2√Σ(1/d²); the earlier continuum Nyquist π√Σ(1/d²) was unnecessarily
restrictive. The RK3 margin is still `field_cfl`×√3. Physical wave and
resistive-decay comparisons check the resulting time integration.

Optional `thermal={...}` selects `EnergyElectrons` by default: conservative
internal-energy transport, adjoint pressure compression, constant/Spitzer η,
and explicit wall/empty-node losses. `closure='entropy'` retains the earlier
experimental entropy scheme only for reproduction. Actual Boris velocities
and all pushed particles supply bulk work velocity. RK3 Hall velocity/work
and Joule heat use stage weights 1/6, 1/6, 2/3. Ions see the ideal/Hall/pressure field;
resistive electric acceleration and ion friction cancel in this mode.
Legacy isothermal/barotropic runs retain their previous path when omitted.
An explicitly requested `thermal={'backend':'metal', ...}` selects the
experimental `MetalEnergyElectrons` subclass. Only flux/wall-rate evaluation
uses GPU float32; conservative integration, compatible work, fields, vacuum
and particles remain CPU64. The default remains CPU; missing GPU/MLX raises.
`gpu_moments=True` additionally selects normalized float32 atomic CIC moments
after the common CPU64 initialization. It is a separately measured prototype:
atomic deposition no longer has the CPU64 gather/deposit adjoint precision.
Particle integration and absorption still use CPU64 positions and velocities.
`thermal={'fused_arrays':True, ...}` optionally fuses CPU64 electron stencils
on either backend. Geometry and the field/RK/particle/circuit steps are not
changed. Edge resistivity is borrowed until the next electric-field call;
stage heating and pressure-power arrays remain owned until the RK sum.
The acoustic bound uses p/max(n,n_floor) for the energy closure, matching its
force denominator; it does not change physical pressure or temperature.

## Geometry
A Cartesian box with three cell types, set by each cell's centre:
- **Gas:** inside the tube, r < r_w and |z| < L. Ions live here.
- An optional `inner_radius` excludes the coaxial inner tube from gas. Its
  interior is vacuum, with no fluid density or pressure; ions entering it
  are absorbed by the same staircase-cell rule. `lost.inner_wall` is a subset
  of radial wall losses. Loss surface classification uses the entered cell's
  centre, matching the staircase geometry, rather than the ideal cylinder.
- **Vacuum gap:** between the tube and the coil, and beyond the tube ends up to
  the conducting end planes. The end planes are at least one cell beyond the
  ends, so the gap closes over the ends.
- **Conductor:** the coil at r ≥ r_c, the end planes, and everything outside
  them. At least two layers line the box, so the periodic indexing of the
  kernels never couples real cells.

**`periodic(device)`** with `device_box(..., periodic_length=L)` gives an
infinitely long tube, periodic in z, for equilibrium and wave tests.

## Fields
**Faraday's law.** Faces touching a gas cell are advanced with the Ohm's-law E
on their edges. E is used only on edges that touch a gas cell and do not lie on
a conductor surface; on a conductor surface tangential E = 0. Conductor faces
keep their normal B: a perfect conductor, the flux conserver of the crowbarred
coil, as in the 2D model.

**Exact vacuum in the gap.** J = 0 there, so B = −∇χ.
- χ solves div B = 0 in every gap cell (a Neumann problem with sparse LU).
- The data come from the gas-side faces and the conductor faces.
- It is applied after every RK3 stage.

An optional `source(t)` gives prescribed, divergence-free coil B on faces.
The vacuum projector acts on B minus this source, so physical coil currents
inside the inner tube are retained. On faces not touching gas, the response
boundary is fixed and the source follows its waveform at each RK stage time.
All faces of gas cells evolve only by Faraday's law using total-field Ohm's
law; the source is not added directly to the plasma. This is a prescribed
external-current boundary, not a self-consistent 3D circuit. No rotational
averaging, mirrored sector or azimuthal filtering is used in evolution.

**Annular vacuum topology.** An annular gas volume is a solid torus; its
connected vacuum has a meridional circulation harmonic. A potential jump
across the inner-vacuum midplane disk, projected to zero divergence and
zero boundary-normal field, supplies this discrete harmonic. Its coefficient
is retained by projection. Gap faces advance by Faraday using Ohm's-law E
on gas edges and prescribed −∂A_source/∂t elsewhere. Thus the circulation
evolves with the boundary EMF instead of being erased by a single-valued
scalar-potential solve. This does not include a rod's Bθ circulation.

**Harmonic field in the periodic tube.** There the gap is a ring around the
plasma and carries one harmonic field: a uniform B_z, the axial flux through
the ring.
- Its gap z-faces are advanced by Faraday's law, with E = 0 on gap-interior
  edges and the Ohm's-law E on the plasma surface.
- The ring flux therefore changes by the line integral of E around the plasma.
- The projection keeps their mean.
- A circulation harmonic (B_θ ∝ 1/r) is not included; it needs an axial current
  through the ring.

## Ions
- **Absorption:** an ion that enters a non-gas cell is absorbed at the
  staircase tube wall or end wall. Its count and energy are recorded (`lost`:
  ends against side wall).
- **Wall-node density:** node density and flux are divided by the gas fraction
  of the node's eight cells. A uniform plasma then has uniform node density up
  to the wall, to 2e-14, like the 2D model's half-area wall nodes. Without this
  the wall nodes read half or a quarter of the density, and the Hall time step
  collapsed (197 instead of 13 field substeps).

`magnetic_energy` integrates only gas and vacuum dual face volumes, with
half weight at the conductor interface. Field values inside conductor cells
serve boundary indexing and are excluded from physical stored energy. Old
runs using whole-box energy have a different constant offset for fixed
boundaries; their time differences are still comparable in that special case.
`response_magnetic_energy` applies the same physical weights to B−B_source;
the coupled-driver load budget uses this response energy and counts every
declared absorbed/cutoff loss. It is separate from total-field stored energy.

## Verification (`examples/hybrid3d_verification.py`)
**Gap projection.** A tilted uniform field (0°, 17°, 69°) is reproduced on all
gap faces to 1e-11 after the gap faces are scrambled; div B in the gap is
≤ 2e-10 (× B/dx).

**Compressional (breathing) eigenmode.** A cold uniform plasma column fills the
tube, with exact vacuum to the conductor. With k_z = 0 the Hall term drops out
and ω = x v_A/a, where

    2 J1(x)/x + ((b² − a²)/a²) J0(x) = 0

a and b are the staircase radii of equal area. Run at 1e22 m⁻³, 0.2 T, with
ω/Ω_i ≈ 0.05, so the radial-only start is the eigenmode:

| Case | dx | ω error | amplitude / exact |
|---|---|---|---|
| Gap to the coil (b = 16.4 cm), matched current support | 10 mm | +6.62% | 0.994 |
| Gap to the coil (b = 16.4 cm), matched current support | 7.5 mm | +3.01% | 0.982 |

These full two-period runs lose no particles and change total energy by
6.72e−8 and 3.99e−8 relative; evidence is in
`outputs/hybrid3d_boundary_mode_check.json`. The former +0.3%/−1.2%
measurements used unmatched wall-node current normalization and are
historical. Staircase interfaces are first order; finer spatial convergence
is still needed before using a 20 mm full shot quantitatively.

**Bugs this test found** (all fixed):
1. Edges on a gas–conductor surface carried Ohm's-law E; the frequency was 41%
   high.
2. With no end-gap cell, the gap ring was capped by conductor and its flux was
   frozen.
3. In a periodic tube the ring's axial flux was dropped, which zeroed the gap
   B_z.
