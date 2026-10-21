# hybrid3d.py

The 3D hybrid core: kinetic full-orbit ions and a massless, quasi-neutral
electron fluid, with Darwin fields on a periodic Yee box. This is milestone M1
of [docs/HYBRID_3D.md](../docs/HYBRID_3D.md).

## Model
**Ions** are macro-particles (weight w = physical ions each).
- Positions sit at integer times, velocities at half steps.
- Boris push with fields gathered trilinearly from nodes.
- Density and flux are deposited by CIC to nodes. The kernels are numba,
  parallel over threads, with per-thread deposit buffers.

**Electrons** follow Ohm's law:

    E = −u_i×B + J×B/(n e) − ∇p_e/(n e) + η J − η_H ∇²J,   J = ∇×B/μ0

- The ideal, Hall and pressure terms are formed at nodes (B, J averaged from
  faces and edges) and averaged to edges.
- η J and the declared hyper-resistivity η_H act on the edge J directly.
- The particles see the node field plus η J.
- p_e is isothermal (`te_ev`) or adiabatic (`gamma_e`).
- `n_floor` is a declared floor on n in Ohm's law.

**Fields.** Faraday's law on faces (Yee: B on faces, E and J on edges), so
∇·B stays at round-off.

## Step
1. E^n from B^n, n^n, and the flux extrapolated 1.5 F^{n−½} − 0.5 F^{n−3/2}.
2. Boris push.
3. Deposit n^{n+1} at x^{n+1}, and the flux at the midpoint positions with
   v^{n+½}.
4. B advanced by SSP-RK3 subcycles, with n linear in time and the flux held at
   F^{n+½}.

The subcycle count comes from the grid-scale bound
ω = k v_A + k² B/(μ0 n e) + (η/μ0) k² + (η_H/μ0) k⁴ (+ k c_s) with
k = π √(Σ 1/d²). It keeps ω dt_B ≤ `field_cfl`·√3, the RK3 imaginary-axis limit.

## Verification (`examples/hybrid3d_verification.py`, 2026-10-20)
**Boris orbit** (Ω dt from 2π/16 to 2π/128, 10 periods):
- The phase error falls 0.79 → 0.20 → 0.050 → 0.013 rad, second order.
- It equals the exact Boris error, 10·2π·(Ω dt)²/12.
- Speed is conserved to 1e-15.

**Parallel waves** in a cold uniform plasma (H⁺, 1e20 m⁻³, 0.1 T), against the
exact cold-ion hybrid dispersion ω² = k²v_A²(1 + ω/Ω_i). Frequency error at
16 / 32 / 64 cells per wavelength:

| k d_i | Whistler | Ion-cyclotron wave |
|---|---|---|
| 0.25 | −5.7% / −1.4% / −0.35% | +4.6% / +1.1% / +0.27% |
| 0.5 | −3.3% / −0.80% / −0.19% | +2.0% / +0.48% / +0.10% |
| 1 | −2.0% / −0.50% / −0.13% | +0.78% / +0.16% / +0.025% |
| 2 | −1.5% / −0.38% / −0.096% | +0.25% / +0.051% / +0.011% |

All of these are second order.
- Total energy (ion kinetic + magnetic; cold electrons, η = 0) drifts by
  ≤ 1.2e-3 of the wave energy over two periods.
- ∇·B ≤ 3e-16 (relative to B/dx).
- The device kernels in `hybrid3d_device.py` reproduce these numpy operators
  to round-off (J identical, E 1e-16, dB/dt 3e-16).
