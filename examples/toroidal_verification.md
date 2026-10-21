# toroidal_verification.py

Verifies the 2.5D (v_θ, B_θ) and Hall options of `ConservativeVacuumMHD`
against exact solutions that are independent of the solver.

```bash
python -m examples.toroidal_verification OUT CASE NR [NR ...] [key=value]
```

`CASE` is `torsional`, `rotor`, `alfvenic`, `whistler` or `rod`. Each writes
`OUT/CASE.json`.

**Setup.** Every case uses the formation settings (staggered energy,
variational boundary Poynting flux, RK3, `ct` magnetic-energy flux, local
trace bounds) without dual energy. The gas cylinder has insulating walls and
ends, and the exact vacuum lies outside. Reference cell averages use 4×4 Gauss
points per cell. Only the torsional wave's starting amplitude is projected.

## Cases
1. **torsional.** A standing torsional Alfvén wave in a uniform B_z.
   - With B_r = 0 every field line is independent, so
     v_θ = V f(r) cos(ks) cos(ωt) and
     B_θ = −√(μ0ρ) V f(r) sin(ks) sin(ωt) are exact for any f with f(R) = 0.
   - s = z + L, k = π/2L, ω = k v_A: free ends, B_θ = 0 at both.
   - f = x(1−x²)³, so the odd reflection at the wall is C².
   - V = 1e-3 v_A and β = 1. Errors are measured after one period.
2. **rotor.** An isothermal rigid rotor, v_θ = W r with
   ρ = ρ0 exp(W²r²/2c²), in a uniform B_z: Mach 0.5 at the wall, 3 sound
   crossings. This tests the centrifugal source against the pressure gradient.
3. **alfvenic.** A nonlinear stationary state with no poloidal field.
   - v = B/√(μ0ρ) along B_θ = B1 f(r) cos³(πz/2L), with
     p = P0 − B_θ²/(2μ0).
   - Hoop stress cancels the centrifugal force and the total pressure is
     uniform.
   - This tests the B_θ magnetic pressure, hoop and centrifugal terms together
     at Mach ~1 and β ~1, over 2 Alfvén crossings.
4. **whistler.** The electron-MHD evolution of a toroidal-field packet in a
   uniform B0, with immobile ions (ion mass 1e6 m_p, so ω_ci/ω ~ 1e-6).
   - **Reference.** The exact linear solution in unbounded space is a
     Hankel–Fourier superposition with ω = D|k||k_z|, D = B0/(μ0 n e):
     b_θ = b0 cos ωt and b_r = −sgn(D)(|k_z|/|k|) b0 sin ωt, evaluated by Gauss
     quadrature. The parallel limit is ω = k² d_i v_A.
   - **Packet.** λ = 4 cm, widths 4 cm (r) and 5 cm (z), amplitude 1e-4 B0. It
     stays e⁻²⁰ below its peak at every boundary for one period.
   - **Measures.** The induced flux at T/4. The B_θ field after one period,
     in the B_θ (energy) norm. The frequency, measured as the scale on D that
     best fits B_θ.
   - `fraction=` sets `hall_whistler_fraction`.
5. **rod.** A current threads an annular chamber (r 4–10 cm, half-length
   30 cm) through `toroidal_boundary` = μ0 I(t)/(2π), with I ramping to 10 kA
   as sin² over 1 µs.
   - B_θ diffuses resistively (η = 2e-3 Ω m) into a heavy, stationary annulus
     (1 kg/m³, 1e5 Pa) from both walls and the ends.
   - **Reference.** At the midplane, a = rB_θ obeys
     ∂a/∂t = (η/μ0) r ∂_r((1/r)∂_r a) with a = μ0 I/(2π) at both walls. It is
     solved independently by second-order finite differences on 4000
     intervals (Radau, rtol 1e-11).
   - **Energy theorem.** The separately counted boundary B_θ Poynting work must
     equal the toroidal magnetic energy plus the Ohmic heat.

## Results (2026-10-20)
See `plasma_bench/conservative_coupled.md` (2.5D section) for the tables.

- **Torsional, rotor, Alfvénic.** All converge at second order or better in L1.
- **Rod** (12, 24, 48 cells across the 6 cm gap). The midplane rB_θ error,
  relative to the final value, is 4.6e-3, 1.2e-3 and 3.2e-4 (second order) at 3 µs. The energy-theorem
  residual over the boundary work is −1.4e-5, −1.2e-6 and −1.3e-7. The plasma moves less than
  0.1 m/s.
- **Whistler.** The frequency error falls at second order. A whistler
  stability scan at λ/7 (fraction 1.0, 1.5, 2.0 over 4 periods) is stable
  through 1.5 and grows at 2.0.
- **Energy.** Fluid-energy residuals are at roundoff. The Hall work is
  ≤ 1e-25 J.

These are the results with the wall-cell Ω slope set to zero
(`mhd_reconstruction.md`). That change leaves every case above unchanged to the
quoted precision.

Runtimes at the finest resolution:
- torsional, 64×128: 285 s;
- rotor, 64×128: 260 s;
- whistler, 144×374: 450 s;
- rod, 48×120: 1130 s.
