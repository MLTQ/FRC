# hybrid3d_arrays.py

Float64 loop fusion and scratch reuse for the existing conservative electron
energy operator. This module changes its evaluation, not its equation,
geometry, pressure work, first-order upwind diffusion or positivity timestep.
Numba uses no fast math. Compression, pressure power, edge eta and the RHS
parallelize independent nodes. The measured compiler lowers the explicitly
unrolled transport and Joule loops serially and emits performance warnings;
those warnings do not indicate a failed calculation. Reported timings include
that behavior. There are no parallel reductions or atomics here.
Parallel induction indices are explicitly signed before periodic neighbor
arithmetic, avoiding unsigned underflow and mixed tuple-index promotion.
Periodic wraps use conditional neighbors rather than integer division in
the inner loop; three-component loops are explicitly unrolled for SIMD.

- `compression`: central divergence of gas-volume velocity, in x/y/z order.
- `hall_power`: the same central pressure gradient and J·E work without
  constructing three gradient/product arrays. Its result is owned, because
  RK stages retain all three work values until their final weighted sum.
- `eta_to_edges`/`joule_heat`: preserve the existing edge average and
  volume-sum-preserving heating transfer. Spitzer evaluation itself remains
  NumPy. Heating results are also owned until their RK combination.
  Currents and masks are passed as component arrays without stacking copies.
- `transport`: the same upwind link fluxes, outgoing coefficient and wall
  enthalpy rate as `EnergyElectrons._transport_state`.
- `energy_rhs`: writes the existing conservative weighted-energy update,
  pressure work and wall outflow into reusable float64 arrays.
- `EnergyArrays`: one workspace per electron closure. Returned arrays are
  borrowed and overwritten by the next corresponding call. No reentrancy.

The base closure retains all positivity checks, installation of updated U,
energy sums, boundary accounting and numerical heat-removal counters.
Metal transport remains an independent provider; its fields are compatible
with the fused CPU64 RHS. No geometric vacuum or force/gather code changes.

| Dependent | Contract |
| --- | --- |
| `EnergyElectrons` | Optional `fused_arrays` workspace, same grid for its lifetime |
| `MetalEnergyElectrons` | Metal supplies fluxes; CPU64 compression/RHS may be fused |

Link/wall geometry is static for the workspace lifetime. Moving a plasma/vacuum
mask requires rebuilding it along with the closure geometry; cached link/wall
arrays are not invalidated automatically. The moving-mask prototype keeps this
option disabled.

Verified CPU64 states match the pre-change implementation exactly at0.75µs
on20×20×40mm and at0.03µs on10mm with1.29million ions. Existing independent
translation/contact/heating comparisons retain their errors and conservation.
Metal whole-step replay uses about19% less time in the measured intervals;
CPU fine uses19–31% less under varying host load. This is not a full-shot
timing or a cure for the physical low-density losses/precision sensitivity.
Evidence and exact source overlays: `outputs/hybrid3d_array_benchmark/`.
Opt-in remains explicit so the running frozen campaign retains its source.
