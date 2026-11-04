# hybrid3d_metal_energy.py

`MetalEnergyElectrons` is an opt-in `EnergyElectrons` subclass. It replaces
only upwind face-flux and wall-rate evaluation with actual Metal float32
kernels. Stored U, conservative divergence, the force-compatible negative
adjoint compression, retained Joule heat and all budget counters remain
float64. There is no total energy, temperature or mode rescale.
The subclass overrides only the transport provider; the shared base class
owns the update and counters, so future closure fixes apply to both backends.
The shared positive-energy transport subcycles use the GPU's fused outgoing
rate estimate, CPU64 wall/expansion combination and the declared 0.8 CFL.

The equations, density cutoff, heat-removal prescription, absorbing geometry
and timestep criteria are those in `hybrid3d_energy_electrons.md`. Face fluxes
and wall rates differ by float32 roundoff. Exact exterior vacuum and Faraday
curl are unchanged. CPU64 sums retain flux cancellation and the pressure-work
contract instead of relying on a long float32 energy accumulation.

Moving-contact, independent constant/Spitzer heating and actual coupled-bank
0–0.75 µs replays pass; corrected-CFL evidence and the exact source archive
are in `outputs/hybrid3d_gpu_cfl_compression/`.
Transport-only endpoint energy differences from CPU are about 1e−9 J in the
physical replay. Full-shot GPU accuracy/performance remain unmeasured.
Numerical low-density losses of the original closure remain physical accuracy
limitations. Each source-frozen queue retains its explicitly selected backend;
the CPU control remains CPU64 while the new six-shot queue uses Metal.
