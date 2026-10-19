# Axis-regular limited reconstruction

`Reconstruction` fits annular averages using x=r²/2. Even variables use a+b x;
the regular odd components vr and Br use r(a+b x), with their exact volume
moments rather than the midpoint radius. Nonuniform-coordinate MC slopes and
periodic axial MC increments limit oscillations. Boundary derivatives use a
limited three-point one-sided estimate; a two-point secant would put the
derivative at the wrong location and reduce wall-force accuracy. Face states
reflect normal velocity/field at the slip wall.
The axis odd trace is exactly zero. This is a derived second-order method,
not a claim to implement every scheme in Mignone (2014).

Total pressure is reconstructed instead of thermal pressure. The geometric
source integrates the same linear-in-x total pressure against dr exactly.
This preserves the special constant-total-pressure, linear-in-x axial-field
equilibrium; it does not guarantee arbitrary FRC equilibrium preservation.
Face/corner thermal pressure is recovered using reconstructed magnetic fields.
By default, scalar density slopes contract toward the positive cell density
when needed. A separately limited thermal-pressure polynomial supplies positive
auxiliary traces. One weight per cell blends the original pressure toward that
polynomial just enough to make every face/corner evaluation pressure positive.
The source integrates the same blended total pressure, including exact Gauss
integration of reconstructed B². Magnetic traces and cell averages are unchanged.
Activations and weights are exposed through `statistics`; `positivity=False`
reproduces the unbounded reference and its explicit failures. This guarantees
positive evaluation traces, not a theorem for the full CT time update. Invalid
states still raise, with no cell floor, energy repair or hidden fallback.

`staggered_corner_fields` reconstructs face-averaged Bz radially and Br axially,
retaining the independent two-state magnetic jumps needed by full LLF UCT.
Cell traces and staggered traces must both use reconstruction; reconstructing
fluid quantities alone leaves the induction update first order.

## Optional local trace bounds (2026-10-19)
`Reconstruction(..., local_bounds=True)` replaces the 1e-10 relative trace
margin for density and thermal pressure with the minimum of the 3×3
neighbouring cell averages, using edge replication at array edges. Slopes
contract about the cell mean, and the pressure blend uses the same bound.

**Why:** with the default margin, a face or corner density can legally fall
to 1e-10 of its cell value. In the NUCTE-III formation such a wall-corner
state (1.1e-16 kg/m³ beside cells at 1.1e-6) gave a fast speed of 6.2e9 m/s.
The LLF UCT dissipation then produced E_θ ≈ 7.8e7 V/m in one stage, and the
run blew up within two steps. The time step, set from cell averages, cannot
see trace speeds.

**Cost:** min-bound limiters clip smooth extrema. On the coil pulse against a
64-cell reference, the relative L1 density error is 1.22% bounded vs 0.74%
unbounded at 16 cells, and 0.487% vs 0.474% at 32. The standing-shock
pressure error changes by 0.15%. The default (`local_bounds=False`) is
unchanged.
