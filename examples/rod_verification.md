# rod_verification.py

End-to-end comparisons for the annular rod-current path. These are analytic
verification fixtures, not FRC formation experiments. No unit-test suite is used.

`python -m examples.rod_verification OUT KIND [nr ...]` writes the case JSON:
- `static`: the 40 kA current-free Bθ=μ0I/(2πr), uniform Bz, density and
  pressure must stay at rest for 5 µs. Measures flow, rBθ, pressure, work and
  energy residuals at successive times and resolutions.
- `static_first`: the same equilibrium through the first-order admissibility
  baseline; necessary because local positivity fallback uses that path.
- `vacuum`: a rod bank drives a highly resistive gas. Compares I(t) with an
  independently integrated plasma-free circuit and Bθ with μ0I(t)/(2πr).
- `closure RUN ...`: computes the stored-run bank/line/gas boundary-work
  budget. Conservation alone does not establish a force-balanced solution.

`case`, `static`, `vacuum` and `closure` retain the formation setup's geometry
and parameters. The static field uses g=K/⟨r²⟩, matching the solver's discrete
current-free energy convention. Changing this file does not validate historical
rod formation runs; those must be rerun after solver corrections.

## Measured correction
The retained pre-fix 12×72 run reached 12.25% rBθ error, 41 km/s poloidal
speed and 67 J of spurious boundary work by 5 µs. After the analytic-vacuum
reconstruction change, 12×72 and 24×144 runs of both orders give:

| Quantity, through 5 µs | Maximum |
|---|---:|
| Relative rBθ error | 6.7e-16 |
| Relative pressure change | 6.7e-15 |
| Poloidal speed | 5.6e-10 m/s |
| Absolute boundary work | 4.4e-14 J |
| Absolute energy residual | 2.7e-11 J |

Evidence: `outputs/rod_verification_corrected/{static,static_first}.json`.
The independent transient annular-diffusion comparison (`toroidal_verification
OUT rod 12 24`) gives field error/K0 0.0046201 → 0.00122484, and relative
energy-theorem residual 1.11e-5 → 1.11e-6. Its 4000/8000-interval reference
differs by at most 3.1e-9 K0. The no-linked-current torsional control retains
magnetic L1 errors 8.25e-4 → 2.03e-4 at 16/32 cells. These verify the solver
correction; the interrupted six-case physical campaign remains incomplete.
