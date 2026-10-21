# hybrid3d_device.py

The 3D hybrid model in the 12" device: milestone M2 of
[docs/HYBRID_3D.md](../docs/HYBRID_3D.md). It shares the particle kernels and
Ohm's law of [hybrid3d.py](hybrid3d.md). The field kernels are numba versions,
identical to the numpy operators to round-off.

## Geometry
A Cartesian box with three cell types, set by each cell's centre:
- **Gas:** inside the tube, r < r_w and |z| < L. Ions live here.
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
| Gap to the coil (b = 16.4 cm) | 10 mm | +0.3% | 0.96 |
| Gap to the coil (b = 16.4 cm) | 7.5 mm | −1.2% | 0.98 |
| Conductor at the wall (b = a) | 10 mm | −4.0% | 0.91 |
| Conductor at the wall (b = a) | 7.5 mm | −3.6% | 0.95 |

The conductor-at-wall case converges slowly: a staircase conductor with
E_t = 0 is first order.

**Bugs this test found** (all fixed):
1. Edges on a gas–conductor surface carried Ohm's-law E; the frequency was 41%
   high.
2. With no end-gap cell, the gap ring was capped by conductor and its flux was
   frozen.
3. In a periodic tube the ring's axial flux was dropped, which zeroed the gap
   B_z.
